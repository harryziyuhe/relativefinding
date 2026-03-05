
import sys
from tqdm import tqdm

# Input files
geno_file = "chr10_LWK_genotypes.txt"
sample_file = "chr10_LWK.samples"
ped_out = "chr10_LWK_haploid_py.ped"

with open(sample_file) as f:
    samples = [line.strip() for line in f if line.strip()]

n_samples = len(samples)

hap0_alleles = {s: [] for s in samples}
hap1_alleles = {s: [] for s in samples}

with open(geno_file) as f:
    for line in tqdm(f):
        cols = line.rstrip('\n').split('\t')
        chrom, pos, snp_id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
        gts = cols[5:]

        if len(gts) != n_samples:
            sys.stderr.write(f"Mismatch: {len(gts)} genotypes but {n_samples} samples\n")
            sys.exit(1)

        for s, gt in zip(samples, gts):
            if gt in ("./.", ".|."):
                hap0_alleles[s].append(("0", "0"))
                hap1_alleles[s].append(("0", "0"))
                continue

            if '|' in gt:
                a, b = gt.split('|')
            elif '/' in gt:
                a, b = gt.split('/')
            else:
                a = b = gt

            def idx_to_base(idx):
                if idx == '.':
                    return "0"
                i = int(idx)
                if i == 0:
                    return ref
                elif i == 1:
                    return alt
                else:
                    return alt

            base_a = idx_to_base(a)
            base_b = idx_to_base(b)

            hap0_alleles[s].append((base_a, base_a))
            hap1_alleles[s].append((base_b, base_b))

with open(ped_out, "w") as out:
    for s in samples:
        fid0 = s
        iid0 = f"{s}_H0"
        fid1 = s
        iid1 = f"{s}_H1"

        prefix0 = [fid0, iid0, "0", "0", "0", "-9"]
        prefix1 = [fid1, iid1, "0", "0", "0", "-9"]

        geno0 = [allele for pair in hap0_alleles[s] for allele in pair]
        geno1 = [allele for pair in hap1_alleles[s] for allele in pair]

        out.write(" ".join(prefix0 + geno0) + "\n")
        out.write(" ".join(prefix1 + geno1) + "\n")

print(f"Wrote haploid-ish PED: {ped_out}")



