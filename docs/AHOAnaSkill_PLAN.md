# AHOAnaSkill — Implementation Plan

> Skill 名称: **AHOAnaSkill**
> 目的: 对不对称氢化 (Asymmetric Hydrogenation, AHO) 实验数据做结构-活性关系 (SAR) 分析
> 形态: Claude Code skill, conda 隔离, agent 主导分析
> Plan 状态: ✅ 已收敛, 待开工

---

## 0. 设计原则 (贯穿全文)

1. **Agent 主导分析** — Python 模块只负责数据管道与"客观计算", 规律总结全部由 agent 完成
2. **基线固定, 分析灵活** — 标准描述符 (RDKit) 作底子; 结构分析走原语 + 临时脚本
3. **累积式扩展** — agent 在每次分析中产出的 ad-hoc 特征持久化, 下次可直接复用
4. **可追溯** — 每个假设都要留痕 (脚本 + 证据图 + 数据 + 结论)
5. **依赖隔离** — conda env `AHOAnaSkill` 与系统/其他项目环境完全隔离
6. **统一入口** — 所有调用走 `aho` 命令, 不允许直接 `python scripts/xxx.py`

---

## 1. 数据规格

### 1.1 CSV (9 列固定顺序)
```
DATA_ID, CAT_NAME, SOL_NAME, PRO_R_NAME, PRO_S_NAME, REA_NAME, TEMP, PRESSURE, EE
```

| 字段 | 类型 | 单位 | 备注 |
|---|---|---|---|
| `DATA_ID` | str | — | 主键, 幂等性靠它 |
| `CAT_NAME` | str | — | 如 `CAT-53`, 对应 `<sdf-dir>/CAT-53.sdf` |
| `SOL_NAME` | str | — | 溶剂分子名 |
| `PRO_R_NAME` | str | — | R 构型产物 |
| `PRO_S_NAME` | str/empty | — | S 构型产物, **可为空** |
| `REA_NAME` | str | — | 反应底物 (待氢化物) |
| `TEMP` | float | **K** | 温度 |
| `PRESSURE` | float | **bar** | 氢气压力 |
| `EE` | float | — | **-1 .. 1 有符号小数**; 正→R 主, 负→S 主 |

### 1.2 SDF 目录
- **扁平结构** — 全部 `.sdf` 文件在同一目录
- 文件名 = `<NAME>.sdf`, 与 CSV 中 `*_NAME` 列严格一致 (大小写敏感)
- SDF 内必须含 `SMILES` 属性字段 (`mol.GetProp("SMILES")`)

---

## 2. 整体架构 (三层)

```
┌─────────────────────────────────────────────────────────────┐
│ Layer 3: Ad-hoc analysis (Agent-authored, per-task)        │
│   • agent 临时写脚本 → 调 L2 原语 → 落 features.family='custom' │
│   • 持久化到 DB; 元数据进 custom_features_meta              │
└─────────────────────────────────────────────────────────────┘
                            ↑ uses
┌─────────────────────────────────────────────────────────────┐
│ Layer 2: AAR_MOL — 结构查询原语 (按需调用, 不预计算)         │
│   • SMARTS / scaffold / MCS / 相似性 / 原子级报告              │
│   • xyzrender 渲染 (含 R/S 自动标注 / 立体叠合 / VdW / 高亮)   │
└─────────────────────────────────────────────────────────────┘
                            ↑ also feeds
┌─────────────────────────────────────────────────────────────┐
│ Layer 1: AAR_FEAT — 基线特征 (固定, 仅 RDKit)                │
│   • rdkit_desc 全集 (~210) + counts + 官能团 SMARTS + Morgan │
│   • 全数据集统一可比的基线                                    │
└─────────────────────────────────────────────────────────────┘
```

**关键: AAR_FEAT 是地板, 不是天花板**. 任何结构性 (位阻 / 邻位 / 配位 / 立体环境) 的判断, agent 必须走 Layer 2 + 临时定义 Layer 3 特征.

---

## 3. 目录结构

```
~/.claude/skills/AHOAnaSkill/
├── SKILL.md                  # agent 入口 (workflow + 调用规约)
├── README.md                 # 人类向快速上手
├── environment.yml           # conda 环境定义
├── install.sh                # 一键安装
├── uninstall.sh              # 干净卸载
├── bin/
│   └── aho                   # 唯一对外命令; conda run 透传
├── scripts/
│   ├── constants.py          # SMARTS 字典 / 金属表 / 列名常量
│   ├── aar_mol.py            # SDF 读取 + StructQuery + Renderer
│   ├── aar_db.py             # SQLite CRUD + import-csv + run-custom-script
│   ├── aar_feat.py           # Layer 1 基线特征
│   ├── ingest.py             # 端到端: CSV+SDF → DB → 基线特征
│   └── check_env.py          # 环境自检
├── data/
│   ├── aho.sqlite            # 持久数据库
│   └── sdf_index.json        # name→path 缓存
└── reports/
    ├── figures/              # xyzrender 渲染输出 (SVG/PNG)
    ├── scripts/              # Layer 3 ad-hoc 脚本归档
    └── YYYYMMDD-HHMM_<topic>.md
```

---

## 4. 数据库 Schema (SQLite)

```sql
CREATE TABLE molecules (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    name          TEXT NOT NULL,
    role          TEXT NOT NULL CHECK(role IN ('REA','CAT','SOL','PRO_R','PRO_S')),
    smiles        TEXT,                -- 优先 SDF 属性, 兜底 MolToSmiles
    inchi_key     TEXT,
    sdf_path      TEXT,
    n_atoms       INTEGER,
    mw            REAL,
    created_at    TEXT,
    UNIQUE(name, role)
);

CREATE TABLE reactions (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    data_id         TEXT UNIQUE NOT NULL,
    reactant_id     INTEGER REFERENCES molecules(id),
    catalyst_id     INTEGER REFERENCES molecules(id),
    solvent_id      INTEGER REFERENCES molecules(id),
    product_r_id    INTEGER REFERENCES molecules(id),
    product_s_id    INTEGER REFERENCES molecules(id),     -- 可空
    temperature_k   REAL,
    pressure_bar    REAL,
    ee              REAL CHECK(ee >= -1 AND ee <= 1),
    source          TEXT,
    notes           TEXT,
    created_at      TEXT
);

CREATE TABLE features (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    molecule_id     INTEGER NOT NULL REFERENCES molecules(id),
    family          TEXT NOT NULL,    -- rdkit_desc / rdkit_count / rdkit_fg /
                                       -- rdkit_morgan / custom_<topic>
    feature_name    TEXT NOT NULL,
    feature_value   REAL,
    computed_at     TEXT,
    UNIQUE(molecule_id, family, feature_name)
);

CREATE TABLE custom_features_meta (
    family_topic    TEXT PRIMARY KEY,    -- 例: custom_ortho_OMe
    description     TEXT,
    hypothesis      TEXT,                -- 对应的结构假设
    script_path     TEXT,                -- 归档脚本路径
    verdict         TEXT,                -- supported/rejected/inconclusive
    n_molecules     INTEGER,
    created_at      TEXT
);

CREATE INDEX idx_reactions_data_id ON reactions(data_id);
CREATE INDEX idx_features_lookup   ON features(molecule_id, family);
```

**去重策略**:
- `molecules`: UNIQUE `(name, role)` — 同名+同角色复用; 同名跨角色独立
- `reactions`: UNIQUE `data_id` — 幂等
- `features`: UNIQUE `(molecule_id, family, feature_name)` — 防重算重写

---

## 5. Conda 环境

### 5.1 `environment.yml`

```yaml
name: AHOAnaSkill
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.12
  - rdkit
  - pandas
  - numpy
  - scipy
  - pip
  - pip:
    - xyzrender
```

> 锁版本到次要号即可; 完全复现时再 `conda env export > environment.lock.yml`.

### 5.2 `install.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="AHOAnaSkill"

# 1) 前置检查
command -v conda >/dev/null 2>&1 || {
    echo "[x] conda 未安装. 请先装 Miniconda/Anaconda."
    exit 1
}

# 2) 创建/更新环境
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[i] env '$ENV_NAME' 已存在, 更新中..."
    conda env update -n "$ENV_NAME" -f "$SKILL_DIR/environment.yml" --prune
else
    echo "[i] 创建 env '$ENV_NAME'..."
    conda env create -f "$SKILL_DIR/environment.yml"
fi

# 3) 可执行权限
chmod +x "$SKILL_DIR/bin/aho"

# 4) 初始化 DB
conda run -n "$ENV_NAME" python "$SKILL_DIR/scripts/aar_db.py" init \
    --db "$SKILL_DIR/data/aho.sqlite"

# 5) 自检
echo "[i] 环境自检..."
conda run -n "$ENV_NAME" python "$SKILL_DIR/scripts/check_env.py"

# 6) PATH 提示
echo
echo "[✓] 安装完成."
echo ""
echo "推荐把 bin 加入 PATH:"
echo "    echo 'export PATH=\"$SKILL_DIR/bin:\$PATH\"' >> ~/.zshrc"
echo "    source ~/.zshrc"
```

### 5.3 `uninstall.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
ENV_NAME="AHOAnaSkill"
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
read -p "卸载 env '$ENV_NAME' 和 skill 数据? [y/N] " ans
[[ "$ans" =~ ^[Yy]$ ]] || exit 0
conda env remove -n "$ENV_NAME" -y || true
echo "[i] env 已删除; skill 目录保留: $SKILL_DIR"
```

### 5.4 `bin/aho` (统一入口)

```bash
#!/usr/bin/env bash
set -euo pipefail
SELF="$(readlink -f "${BASH_SOURCE[0]}")"
SKILL_DIR="$(cd "$(dirname "$SELF")/.." && pwd)"
SCRIPTS="$SKILL_DIR/scripts"
ENV_NAME="AHOAnaSkill"
_py() { conda run --no-capture-output -n "$ENV_NAME" python "$@"; }
CMD="${1:-help}"; shift || true
case "$CMD" in
    ingest)         _py "$SCRIPTS/ingest.py" "$@" ;;
    db)             _py "$SCRIPTS/aar_db.py" "$@" ;;
    mol)            _py "$SCRIPTS/aar_mol.py" "$@" ;;
    feat)           _py "$SCRIPTS/aar_feat.py" "$@" ;;
    render)         _py "$SCRIPTS/aar_mol.py" render "$@" ;;
    overlay)        _py "$SCRIPTS/aar_mol.py" overlay "$@" ;;
    smarts)         _py "$SCRIPTS/aar_mol.py" smarts-highlight "$@" ;;
    run-custom)     _py "$SCRIPTS/aar_db.py" run-custom-script "$@" ;;
    stats)          _py "$SCRIPTS/aar_db.py" stats ;;
    check-env)      _py "$SCRIPTS/check_env.py" ;;
    help|--help|-h) sed -n '/^### aho commands/,/^###/p' "$SKILL_DIR/README.md" ;;
    *) echo "[x] 未知命令: $CMD; 用 'aho help' 看完整列表"; exit 1 ;;
esac
```

### 5.5 `scripts/check_env.py`

```python
"""Env health check. Exit non-zero if any required package missing."""
import sys, sqlite3

def check(name, fn):
    try:
        v = fn(); print(f"  [✓] {name}: {v}"); return True
    except Exception as e:
        print(f"  [x] {name}: {e}"); return False

print(f"Python: {sys.version.split()[0]}")
ok = all([
    check("rdkit",     lambda: __import__("rdkit").__version__),
    check("xyzrender", lambda: __import__("xyzrender").__version__),
    check("pandas",    lambda: __import__("pandas").__version__),
    check("numpy",     lambda: __import__("numpy").__version__),
    check("scipy",     lambda: __import__("scipy").__version__),
    check("sqlite3",   lambda: sqlite3.sqlite_version),
])
sys.exit(0 if ok else 1)
```

---

## 6. 三个模块 API (设计契约)

### 6.1 `aar_mol.py` — Layer 2 结构原语

```python
class MolReader:
    def read_sdf(path) -> dict:
        # → {name, smiles, inchi, inchi_key, n_atoms, mw, sdf_path}
        # SMILES 先取 SDF 属性, 兜底 MolToSmiles
    def build_index(sdf_dir) -> dict[name, path]
    def validate(name, index) -> Optional[path]

class StructQuery:
    # 子结构 / 官能团
    smarts_count(smiles, smarts) -> int
    smarts_match(smiles, smarts) -> list[tuple]
    has_substructure(smiles, smarts) -> bool
    functional_groups(smiles) -> dict       # 内置 ~40 个 FG SMARTS

    # 骨架 / 相似
    bemis_murcko_scaffold(smiles) -> str
    generic_skeleton(smiles) -> str
    mcs(smi_a, smi_b) -> str
    tanimoto(smi_a, smi_b) -> float
    cluster_by_scaffold(smi_list) -> dict

    # 原子级
    atom_report(smiles) -> list[dict]       # idx/sym/hyb/charge/ring/aromatic/chir/nbr
    stereo_centers(smiles) -> list[dict]
    coordination_env(cat_smiles) -> dict    # 催化剂: metal/donors/ring_sizes/bite_est

    # 量化代理 (纯 RDKit)
    count_ring_systems(smi) -> int
    count_aromatic_atoms(smi) -> int
    longest_aliphatic_chain(smi) -> int
    sterimol_proxy(smi, attach_idx) -> dict  # 简化版 L/B1/B5

class Renderer:  # xyzrender 封装
    render(sdf, out, preset='default') -> path
    render_stereo(sdf, out, label_stereo=True) -> path
    render_pair_overlay(sdf_r, sdf_s, out, colors=('blue','red')) -> path
    render_with_vdw(sdf, out, atoms_region=None) -> path
    render_with_annotations(sdf, out, distances=None, angles=None) -> path
    render_with_smarts_highlight(sdf, smiles, smarts, out) -> path
    render_grid(sdf_list, labels, out, n_cols=4) -> path   # SVG 拼接
    render_rotation_gif(sdf, out) -> path
```

**CLI**:
```
aho mol read --sdf X.sdf
aho mol smarts-count --smi "..." --pattern "[c]-[OH]"
aho mol atom-report --smi "..." --json
aho render --sdf X.sdf --out X.png --preset paton
aho overlay --sdf-a A.sdf --sdf-b B.sdf --out pair.png
aho smarts --sdf X.sdf --pattern "c-OMe" --out hl.png
```

### 6.2 `aar_db.py` — Layer 1 持久层

```python
class DB:
    def init(db_path)
    def upsert_molecule(name, role, sdf_path) -> mol_id
    def upsert_reaction(row: dict)             # 处理一行 CSV
    def import_csv(csv_path, sdf_dir, db_path) -> summary
    def query_reactions(filter: dict) -> list[dict]
    def query_with_features(filter, families) -> wide_table   # 给 agent 的宽表
    def stats() -> dict                          # 数量/分布/EE 直方
    def run_custom_script(script_path)           # 执行 ad-hoc + 落 meta
    def register_custom_meta(family_topic, description, hypothesis, script_path, verdict)
```

**CLI**:
```
aho db init
aho db import-csv --csv X.csv --sdf-dir Y/
aho db query --catalyst CAT-53 --json
aho db stats
aho run-custom --script reports/scripts/X.py
```

### 6.3 `aar_feat.py` — Layer 1 基线特征 (仅 RDKit)

```python
class RDKitFeaturizer:
    def compute_for_molecule(mol_id) -> int   # 入库特征数
    def compute_all_new(family='rdkit_desc')  # 跳过已算

    # 内部:
    def _standard_desc(mol)  -> dict          # Descriptors._descList 全集 (~210)
    def _counts(mol)         -> dict          # n_atoms/n_rings/n_aromatic/n_stereo/fsp3
    def _fg_inventory(mol)   -> dict          # ~40 个官能团 SMARTS 计数
    def _morgan_fp(mol, r=2, n=1024) -> dict  # active bits → sparse 存
```

**CLI**:
```
aho feat compute --family rdkit_desc
aho feat compute --family rdkit_morgan --role CAT
aho feat stats
```

**不含**: xTB / SOAP / ACSF / dscribe / 任何外部二进制依赖.

---

## 7. SKILL.md — agent 工作流程

### 7.1 触发词
不对称氢化 / AHO / AAR / 氢化数据 / 构效关系 / SAR / 这批数据看一下规律

### 7.2 执行步骤

1. **环境检查**: `aho check-env && aho stats`
2. **摄入**: `aho ingest --csv X --sdf-dir Y`
   - 缺 SDF 行**不静默跳过**, 明确告知用户
3. **基线特征**: `aho feat compute --family rdkit_desc`
4. **数据汇出**: `aho db query-with-features --json` → 拿宽表
5. **看结构 (强制)**:
   - 每个 role 至少渲一个代表样本: `aho render --sdf ...`
   - PRO_R/PRO_S 至少一对: `aho overlay --sdf-a ... --sdf-b ...`
6. **提结构假设** (agent 主导, 非固化):
   - 同 REA + 同 CAT, 变 TEMP/PRESSURE/SOL → 条件效应
   - 同 REA, 变 CAT → 催化剂筛选
   - 同 CAT, 变 REA → 底物 scope
   - EE 符号反转 → 专项讨论
   - 邻位/位阻类: 必须先 `aho smarts --pattern ...` 在图上确认子结构位置
7. **量化验证 (Layer 3)**:
   - 写 ad-hoc 脚本: 用 `StructQuery.smarts_count` 等原语跑全表
   - 归档到 `reports/scripts/<ts>_<topic>.py`
   - `aho run-custom --script ...` 执行 → 写入 `features.family='custom_<topic>'`
   - 注册 `custom_features_meta`: hypothesis + verdict (supported/rejected/inconclusive)
8. **统计**: 相关系数 / 分组均值差 / Mann-Whitney U / 异常点 (>2σ)
9. **报告**: `reports/<ts>_<topic>.md`
   - 数据概览表 (n_reactions / 各 role 唯一数 / EE 直方文字版)
   - 每个结论: n + 证据 + 反例 + 引用证据图
   - **Reasoning Trace** 章节: 列出本次产生的所有 custom features, 各自 verdict 与脚本路径
   - 建议的下一步实验

### 7.3 调用规约

```
✅ 必须: aho <subcommand> ...
❌ 不要: python /full/path/scripts/xxx.py ...   (用错 Python / 缺包)

ad-hoc 脚本里 (在 env 内执行) 可以直接:
    from aar_mol import StructQuery, Renderer
    from aar_db import DB
```

### 7.4 渲图准则

| 何时**必须**渲图 | 何时**不要**渲图 |
|---|---|
| 首次见数据集时, 各 role 代表样本 | 已渲过同一对照的下游分析 |
| 任何位阻/邻位/构象类假设前 | 纯计数/比例统计 |
| PRO_R vs PRO_S 任何讨论 | 已经在 features 表里能直接回答的 |
| 报告里每个核心假设 | — |

图存路径: `reports/figures/<ts>_<topic>_<idx>.svg`

---

## 8. 报告模板要点

```markdown
# AHO Analysis Report — <Topic>
- Generated: <ISO timestamp>
- Data scope: <n_reactions> reactions over <n_REA> reactants × <n_CAT> catalysts
- Feature families used: rdkit_desc, custom_<...>, ...

## Dataset overview
| role  | unique count |
| REA   | …            |
| CAT   | …            |
| ...
EE distribution: <文字直方或散点描述>

## Key findings
### Finding 1: <一句话结论>
- Evidence: n=<x>, mean |ee| = ... (group A) vs ... (group B), p=...
- Counterexamples: <count> rows; <names>
- Figure: ![](figures/<...>.svg)

### Finding 2: ...

## Reasoning Trace
| family_topic | hypothesis | verdict | n_mol | script |
| custom_ortho_OMe | 邻位 OMe 降 ee | supported | 23 | reports/scripts/... |
| ...

## Suggested next experiments
- ...
- ...

## Limitations
- <小样本/缺失数据/未验证假设>
```

---

## 9. 实施 Phase 路线

| Phase | 内容 | 交付物 | 验收 |
|---|---|---|---|
| **P0 — Skeleton** | 目录, environment.yml, install.sh, bin/aho, check_env.py, SKILL.md 草稿 | 可装可卸 | `aho check-env` 全 ✓ |
| **P1 — DB + Ingest** | schema, aar_mol.read_sdf, aar_db.import_csv, ingest.py | 端到端入库 | 小样本 CSV+SDF, `aho stats` 数字正确 |
| **P2 — Renderer** | aar_mol.Renderer 全套 xyzrender 封装 | 7 个渲染子命令可用 | 任选 5 个分子出图无报错 |
| **P3 — Baseline feat** | aar_feat.py: rdkit_desc + counts + fg + morgan | features 表满 | 覆盖率 100% |
| **P4 — Agent workflow** | SKILL.md 完整, aho run-custom 链路, custom_features_meta | 端到端 SAR 案例 | agent 自主走完 7.2 全流程 |
| **P5 — Report** | 模板 + Reasoning Trace + 图引用规范 | 可读报告 | 人评通过 |

---

## 10. 待确认事项 (实施前)

1. **测试数据**: 小样本 CSV (10–30 行) + 对应 SDF 文件夹位置?
   - 没有的话, P1 阶段我用合成数据 (5 个简单底物 + 3 个虚拟催化剂) 先打骨架
2. **首次发布位置**: 直接落到 `~/.claude/skills/AHOAnaSkill/`, 还是先在某个 staging 目录?
3. **报告语言**: 中文 / 英文 / 双语?

---

## 11. 参考代码

- **AAReact** (https://github.com/iawnix/AAReact) — featurizer 类、SDF 读取惯例、name→role 约定
- **xyzrender** (https://github.com/aligfellow/xyzrender) — 渲染层全部能力
- **xyz2svg** (AAReact 上游) — xyzrender 内部 SVG 引擎

---

## 12. 状态

✅ Plan 收敛, 设计完整, 准备进入 P0.

下次会话直接 "开始 P0" 即可启动实施.
