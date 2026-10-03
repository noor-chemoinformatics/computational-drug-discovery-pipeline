"""
Calculate enrichment factor for deep docking results.

Usage:
    python calculate_enrichment.py <library.smi> <final_smiles.smi> <actives.smi> <output.txt>
"""
import sys


def calculate_enrichment(library_smi, final_smiles_smi, actives_smi, output_txt):
    with open(library_smi) as f:
        library = [line.strip().split()[0] for line in f if line.strip()]
    with open(final_smiles_smi) as f:
        final = [line.strip() for line in f if line.strip()]
    with open(actives_smi) as f:
        actives = set(line.strip().split()[0] for line in f if line.strip())

    total_compounds = len(library)
    total_actives = sum(1 for smi in library if smi in actives)
    actives_in_final = sum(1 for smi in final if smi in actives)

    fraction_actives_in_final = actives_in_final / len(final) if final else 0
    fraction_actives_total = total_actives / total_compounds if total_compounds else 0
    ef = fraction_actives_in_final / fraction_actives_total if fraction_actives_total else 0
    recall = actives_in_final / total_actives if total_actives else 0

    lines = [
        f"Total compounds in library: {total_compounds}",
        f"Total actives in library: {total_actives}",
        f"Final shortlist size: {len(final)}",
        f"Actives recovered: {actives_in_final}",
        f"Recall: {recall:.1%}",
        f"Enrichment factor: {ef:.2f}",
    ]
    output = "\n".join(lines)
    print(output)
    with open(output_txt, 'w') as f:
        f.write(output + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python calculate_enrichment.py <library.smi> <final_smiles.smi> <actives.smi> <output.txt>")
        sys.exit(1)
    calculate_enrichment(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
