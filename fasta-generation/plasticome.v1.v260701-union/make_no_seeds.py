import csv

# 1. collect seed accessions (skip the empty rows)
seedset = set()
for r in csv.DictReader(open('pazy-seed-subset/putative-curated-seed-set.csv', newline='')):
    if (r.get('aa_sequence') or '').strip() and (r.get('accession') or '').strip():
        seedset.add(r['accession'].strip())

# 2. read master FASTA
def read_fasta(p):
    hdr, buf = None, []
    for line in open(p):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if hdr is not None: yield hdr, ''.join(buf)
            hdr, buf = line[1:], []
        else: buf.append(line)
    if hdr is not None: yield hdr, ''.join(buf)

def header_accs(hdr):                       # 2nd pipe-field, ';'-split
    parts = hdr.split('|')
    return {t.strip() for t in (parts[1] if len(parts) > 1 else '').split(';') if t.strip()}

def wrap(s, w=60):
    return '\n'.join(s[i:i+w] for i in range(0, len(s), w))

# 3. write complement (records whose accessions don't intersect the seed set)
kept = dropped = 0
matched = set()
with open('pazy-seed-subset/plasticome.v1.v260701-union.no-seeds.fasta', 'w') as o:
    for hdr, seq in read_fasta('plasticome.v1.v260701-union.fasta'):
        hit = header_accs(hdr) & seedset
        if not hit:
            o.write(f'>{hdr}\n{wrap(seq)}\n')
            kept += 1
        else:
            dropped += 1
            matched |= hit

print(f'seeds={len(seedset)} kept={kept} dropped={dropped}')
print(f'seed accessions with no match: {sorted(seedset - matched)}')
