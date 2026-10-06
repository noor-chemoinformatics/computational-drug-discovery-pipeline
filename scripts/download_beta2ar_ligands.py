"""
Download known beta-2 adrenergic receptor ligands from ChEMBL.

Usage:
    python download_beta2ar_ligands.py <output.smi>
"""
import sys
from chembl_webresource_client.new_client import new_client


def download_ligands(output_smi):
    target = new_client.target
    target_query = target.search('beta-2 adrenergic receptor')
    targets = [t for t in target_query if t.get('organism') == 'Homo sapiens']

    if not targets:
        print("No human beta2AR target found")
        return

    target_id = targets[0]['target_chembl_id']
    print(f"Using target: {target_id}")

    activity = new_client.activity
    res = activity.filter(target_chembl_id=target_id).filter(
        standard_type__in=['IC50', 'EC50', 'Ki']
    ).filter(standard_value__isnull=False)

    seen = set()
    count = 0
    with open(output_smi, 'w') as f:
        for act in res:
            smi = act.get('canonical_smiles')
            if smi and smi not in seen:
                seen.add(smi)
                f.write(smi + "\n")
                count += 1
                if count >= 200:
                    break

    print(f"Wrote {count} unique SMILES to {output_smi}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python download_beta2ar_ligands.py <output.smi>")
        sys.exit(1)
    download_ligands(sys.argv[1])
