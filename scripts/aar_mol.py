#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Draw, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

from constants import FUNCTIONAL_GROUP_SMARTS, METAL_SYMBOLS, resolve_smarts


class MolReader:
    @staticmethod
    def build_index(sdf_dir: str | Path) -> dict[str, str]:
        sdf_root = Path(sdf_dir)
        if not sdf_root.is_dir():
            raise FileNotFoundError(f"SDF directory not found: {sdf_root}")
        return {path.stem: str(path.resolve()) for path in sorted(sdf_root.glob("*.sdf"))}

    @staticmethod
    def validate(name: str, index: dict[str, str]) -> str | None:
        if not name:
            return None
        return index.get(name)

    @staticmethod
    def load_sdf(path: str | Path) -> Chem.Mol:
        sdf_path = Path(path)
        if not sdf_path.exists():
            raise FileNotFoundError(f"SDF file not found: {sdf_path}")
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
        if mol is None:
            raise ValueError(f"No readable molecule in SDF: {sdf_path}")
        return mol

    @staticmethod
    def read_sdf(path: str | Path) -> dict:
        sdf_path = Path(path).resolve()
        mol = MolReader.load_sdf(sdf_path)
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
        smiles = mol.GetProp("SMILES").strip() if mol.HasProp("SMILES") else ""
        if not smiles:
            smiles = Chem.MolToSmiles(mol_no_h)
        inchi_key = None
        try:
            from rdkit.Chem import inchi

            inchi_key = inchi.MolToInchiKey(mol_no_h)
        except Exception:
            inchi_key = None
        return {
            "name": sdf_path.stem,
            "smiles": smiles,
            "inchi_key": inchi_key,
            "sdf_path": str(sdf_path),
            "n_atoms": int(mol_no_h.GetNumAtoms()),
            "mw": float(Descriptors.MolWt(mol_no_h)),
        }


class StructQuery:
    @staticmethod
    def _mol(smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return mol

    @staticmethod
    def _pattern(smarts: str) -> Chem.Mol:
        pattern = Chem.MolFromSmarts(resolve_smarts(smarts))
        if pattern is None:
            raise ValueError(f"Invalid SMARTS pattern: {smarts}")
        return pattern

    @classmethod
    def smarts_count(cls, smiles: str, smarts: str) -> int:
        mol = cls._mol(smiles)
        return len(mol.GetSubstructMatches(cls._pattern(smarts), uniquify=True))

    @classmethod
    def smarts_match(cls, smiles: str, smarts: str) -> list[tuple[int, ...]]:
        mol = cls._mol(smiles)
        return [tuple(match) for match in mol.GetSubstructMatches(cls._pattern(smarts))]

    @classmethod
    def has_substructure(cls, smiles: str, smarts: str) -> bool:
        return cls._mol(smiles).HasSubstructMatch(cls._pattern(smarts))

    @classmethod
    def functional_groups(cls, smiles: str) -> dict[str, int]:
        return {
            name: cls.smarts_count(smiles, smarts)
            for name, smarts in FUNCTIONAL_GROUP_SMARTS.items()
        }

    @classmethod
    def bemis_murcko_scaffold(cls, smiles: str) -> str:
        mol = cls._mol(smiles)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold) if scaffold is not None else ""

    @classmethod
    def generic_skeleton(cls, smiles: str) -> str:
        mol = cls._mol(smiles)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold is None:
            return ""
        generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return Chem.MolToSmiles(generic) if generic is not None else ""

    @classmethod
    def mcs(cls, smi_a: str, smi_b: str) -> str:
        mols = [cls._mol(smi_a), cls._mol(smi_b)]
        result = rdFMCS.FindMCS(mols, timeout=10)
        return result.smartsString

    @classmethod
    def tanimoto(cls, smi_a: str, smi_b: str) -> float:
        fp_a = AllChem.GetMorganFingerprintAsBitVect(cls._mol(smi_a), radius=2, nBits=2048)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(cls._mol(smi_b), radius=2, nBits=2048)
        return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))

    @classmethod
    def cluster_by_scaffold(cls, smi_list: Iterable[str]) -> dict[str, list[int]]:
        clusters: dict[str, list[int]] = defaultdict(list)
        for idx, smiles in enumerate(smi_list):
            clusters[cls.bemis_murcko_scaffold(smiles)].append(idx)
        return dict(clusters)

    @classmethod
    def atom_report(cls, smiles: str) -> list[dict]:
        mol = cls._mol(smiles)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        report = []
        for atom in mol.GetAtoms():
            report.append(
                {
                    "idx": atom.GetIdx(),
                    "symbol": atom.GetSymbol(),
                    "hybridization": str(atom.GetHybridization()),
                    "formal_charge": atom.GetFormalCharge(),
                    "ring": atom.IsInRing(),
                    "aromatic": atom.GetIsAromatic(),
                    "chirality": atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "",
                    "neighbors": [nbr.GetIdx() for nbr in atom.GetNeighbors()],
                }
            )
        return report

    @classmethod
    def stereo_centers(cls, smiles: str) -> list[dict]:
        mol = cls._mol(smiles)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
        return [{"idx": int(idx), "label": str(label)} for idx, label in centers]

    @classmethod
    def coordination_env(cls, cat_smiles: str) -> dict:
        mol = cls._mol(cat_smiles)
        donors = []
        metals = []
        donor_symbols = {"N", "O", "P", "S", "C"}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in METAL_SYMBOLS:
                metals.append({"idx": atom.GetIdx(), "symbol": symbol})
            if symbol in donor_symbols and atom.GetFormalCharge() <= 1:
                donors.append(
                    {
                        "idx": atom.GetIdx(),
                        "symbol": symbol,
                        "aromatic": atom.GetIsAromatic(),
                        "neighbors": [nbr.GetIdx() for nbr in atom.GetNeighbors()],
                    }
                )
        ring_sizes = sorted({len(ring) for ring in mol.GetRingInfo().AtomRings()})
        return {
            "metals": metals,
            "donors": donors,
            "ring_sizes": ring_sizes,
            "bite_est": len(donors),
        }

    @classmethod
    def count_ring_systems(cls, smiles: str) -> int:
        mol = cls._mol(smiles)
        rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
        systems: list[set[int]] = []
        for ring in rings:
            for system in systems:
                if ring & system:
                    system |= ring
                    break
            else:
                systems.append(set(ring))
        return len(systems)

    @classmethod
    def count_aromatic_atoms(cls, smiles: str) -> int:
        return sum(1 for atom in cls._mol(smiles).GetAtoms() if atom.GetIsAromatic())

    @classmethod
    def longest_aliphatic_chain(cls, smiles: str) -> int:
        mol = cls._mol(smiles)
        carbon_idxs = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetSymbol() == "C" and not atom.GetIsAromatic() and not atom.IsInRing()
        ]
        if not carbon_idxs:
            return 0
        graph = {idx: [] for idx in carbon_idxs}
        carbon_set = set(carbon_idxs)
        for idx in carbon_idxs:
            atom = mol.GetAtomWithIdx(idx)
            graph[idx] = [
                nbr.GetIdx()
                for nbr in atom.GetNeighbors()
                if nbr.GetIdx() in carbon_set
            ]

        def longest_from(start: int, seen: set[int]) -> int:
            lengths = [longest_from(nbr, seen | {nbr}) for nbr in graph[start] if nbr not in seen]
            return 1 + (max(lengths) if lengths else 0)

        return max(longest_from(idx, {idx}) for idx in carbon_idxs)

    @classmethod
    def sterimol_proxy(cls, smiles: str, attach_idx: int) -> dict:
        mol = Chem.AddHs(cls._mol(smiles))
        if attach_idx < 0 or attach_idx >= mol.GetNumAtoms():
            raise ValueError(f"attach_idx out of range: {attach_idx}")
        if AllChem.EmbedMolecule(mol, randomSeed=17) != 0:
            raise ValueError("Could not embed molecule for sterimol proxy")
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        conf = mol.GetConformer()
        origin = conf.GetAtomPosition(attach_idx)
        distances = []
        for atom in mol.GetAtoms():
            if atom.GetIdx() == attach_idx:
                continue
            pos = conf.GetAtomPosition(atom.GetIdx())
            distances.append(math.dist((origin.x, origin.y, origin.z), (pos.x, pos.y, pos.z)))
        if not distances:
            return {"L": 0.0, "B1": 0.0, "B5": 0.0}
        return {"L": max(distances), "B1": min(distances), "B5": max(distances)}


class Renderer:
    @staticmethod
    def _mol_for_draw(path: str | Path) -> Chem.Mol:
        mol = MolReader.load_sdf(path)
        Chem.rdDepictor.Compute2DCoords(mol)
        return mol

    @staticmethod
    def _write_mol(
        mol: Chem.Mol,
        out: str | Path,
        legend: str = "",
        highlight_atoms: Iterable[int] | None = None,
        size: tuple[int, int] = (700, 520),
    ) -> str:
        out_path = Path(out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        highlight = list(highlight_atoms or [])
        if out_path.suffix.lower() == ".svg":
            drawer = Draw.MolDraw2DSVG(*size)
            drawer.drawOptions().addStereoAnnotation = True
            drawer.DrawMolecule(mol, legend=legend, highlightAtoms=highlight)
            drawer.FinishDrawing()
            out_path.write_text(drawer.GetDrawingText(), encoding="utf-8")
        else:
            Draw.MolToFile(mol, str(out_path), size=size, legend=legend, highlightAtoms=highlight)
        return str(out_path)

    @classmethod
    def render(cls, sdf: str | Path, out: str | Path, preset: str = "default") -> str:
        mol = cls._mol_for_draw(sdf)
        return cls._write_mol(mol, out, legend=f"{Path(sdf).stem} ({preset})")

    @classmethod
    def render_stereo(cls, sdf: str | Path, out: str | Path, label_stereo: bool = True) -> str:
        mol = cls._mol_for_draw(sdf)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        return cls._write_mol(mol, out, legend=f"{Path(sdf).stem} stereo")

    @classmethod
    def render_pair_overlay(
        cls,
        sdf_r: str | Path,
        sdf_s: str | Path,
        out: str | Path,
        colors: tuple[str, str] = ("blue", "red"),
    ) -> str:
        mols = [cls._mol_for_draw(sdf_r), cls._mol_for_draw(sdf_s)]
        legends = [f"{Path(sdf_r).stem} ({colors[0]})", f"{Path(sdf_s).stem} ({colors[1]})"]
        out_path = Path(out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        use_svg = out_path.suffix.lower() == ".svg"
        image = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(450, 360), legends=legends, useSVG=use_svg)
        if use_svg:
            out_path.write_text(image, encoding="utf-8")
        else:
            image.save(out_path)
        return str(out_path)

    @classmethod
    def render_with_vdw(cls, sdf: str | Path, out: str | Path, atoms_region: Iterable[int] | None = None) -> str:
        mol = cls._mol_for_draw(sdf)
        return cls._write_mol(mol, out, legend=f"{Path(sdf).stem} VdW proxy", highlight_atoms=atoms_region)

    @classmethod
    def render_with_annotations(
        cls,
        sdf: str | Path,
        out: str | Path,
        distances: list | None = None,
        angles: list | None = None,
    ) -> str:
        mol = cls._mol_for_draw(sdf)
        suffix = []
        if distances:
            suffix.append(f"{len(distances)} distances")
        if angles:
            suffix.append(f"{len(angles)} angles")
        legend = f"{Path(sdf).stem} {'; '.join(suffix)}".strip()
        return cls._write_mol(mol, out, legend=legend)

    @classmethod
    def render_with_smarts_highlight(
        cls,
        sdf: str | Path,
        smiles: str | None,
        smarts: str,
        out: str | Path,
    ) -> str:
        mol = cls._mol_for_draw(sdf)
        query_mol = StructQuery._mol(smiles) if smiles else Chem.RemoveHs(mol)
        matches = query_mol.GetSubstructMatches(StructQuery._pattern(smarts))
        atoms = sorted({idx for match in matches for idx in match})
        return cls._write_mol(mol, out, legend=f"{Path(sdf).stem}: {smarts}", highlight_atoms=atoms)

    @classmethod
    def render_grid(cls, sdf_list: Iterable[str | Path], labels: list[str], out: str | Path, n_cols: int = 4) -> str:
        mols = [cls._mol_for_draw(sdf) for sdf in sdf_list]
        out_path = Path(out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        use_svg = out_path.suffix.lower() == ".svg"
        image = Draw.MolsToGridImage(mols, molsPerRow=n_cols, subImgSize=(360, 300), legends=labels, useSVG=use_svg)
        if use_svg:
            out_path.write_text(image, encoding="utf-8")
        else:
            image.save(out_path)
        return str(out_path)

    @classmethod
    def render_rotation_gif(cls, sdf: str | Path, out: str | Path) -> str:
        mol = cls._mol_for_draw(sdf)
        out_path = Path(out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        image = Draw.MolToImage(mol, size=(520, 420), legend=f"{Path(sdf).stem} static rotation proxy")
        image.save(out_path, save_all=True, append_images=[image] * 3, duration=180, loop=0)
        return str(out_path)


def _json_print(obj: object) -> None:
    print(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AHO molecule utilities")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("read")
    p.add_argument("--sdf", required=True)

    p = sub.add_parser("smarts-count")
    p.add_argument("--smi", required=True)
    p.add_argument("--pattern", required=True)

    p = sub.add_parser("smarts-match")
    p.add_argument("--smi", required=True)
    p.add_argument("--pattern", required=True)

    p = sub.add_parser("has-substructure")
    p.add_argument("--smi", required=True)
    p.add_argument("--pattern", required=True)

    p = sub.add_parser("functional-groups")
    p.add_argument("--smi", required=True)

    p = sub.add_parser("atom-report")
    p.add_argument("--smi", required=True)
    p.add_argument("--json", action="store_true")

    p = sub.add_parser("stereo-centers")
    p.add_argument("--smi", required=True)

    p = sub.add_parser("scaffold")
    p.add_argument("--smi", required=True)
    p.add_argument("--generic", action="store_true")

    p = sub.add_parser("mcs")
    p.add_argument("--smi-a", required=True)
    p.add_argument("--smi-b", required=True)

    p = sub.add_parser("tanimoto")
    p.add_argument("--smi-a", required=True)
    p.add_argument("--smi-b", required=True)

    p = sub.add_parser("coord-env")
    p.add_argument("--smi", required=True)

    p = sub.add_parser("render")
    p.add_argument("--sdf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--preset", default="default")

    p = sub.add_parser("overlay")
    p.add_argument("--sdf-a", required=True)
    p.add_argument("--sdf-b", required=True)
    p.add_argument("--out", required=True)

    p = sub.add_parser("smarts-highlight")
    p.add_argument("--sdf", required=True)
    p.add_argument("--smiles")
    p.add_argument("--pattern", required=True)
    p.add_argument("--out", required=True)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "read":
        _json_print(MolReader.read_sdf(args.sdf))
    elif args.cmd == "smarts-count":
        print(StructQuery.smarts_count(args.smi, args.pattern))
    elif args.cmd == "smarts-match":
        _json_print(StructQuery.smarts_match(args.smi, args.pattern))
    elif args.cmd == "has-substructure":
        print("true" if StructQuery.has_substructure(args.smi, args.pattern) else "false")
    elif args.cmd == "functional-groups":
        _json_print(StructQuery.functional_groups(args.smi))
    elif args.cmd == "atom-report":
        report = StructQuery.atom_report(args.smi)
        _json_print(report) if args.json else print(report)
    elif args.cmd == "stereo-centers":
        _json_print(StructQuery.stereo_centers(args.smi))
    elif args.cmd == "scaffold":
        print(StructQuery.generic_skeleton(args.smi) if args.generic else StructQuery.bemis_murcko_scaffold(args.smi))
    elif args.cmd == "mcs":
        print(StructQuery.mcs(args.smi_a, args.smi_b))
    elif args.cmd == "tanimoto":
        print(f"{StructQuery.tanimoto(args.smi_a, args.smi_b):.6f}")
    elif args.cmd == "coord-env":
        _json_print(StructQuery.coordination_env(args.smi))
    elif args.cmd == "render":
        print(Renderer.render(args.sdf, args.out, preset=args.preset))
    elif args.cmd == "overlay":
        print(Renderer.render_pair_overlay(args.sdf_a, args.sdf_b, args.out))
    elif args.cmd == "smarts-highlight":
        print(Renderer.render_with_smarts_highlight(args.sdf, args.smiles, args.pattern, args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
