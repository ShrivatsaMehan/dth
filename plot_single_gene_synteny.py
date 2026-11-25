#!/usr/bin/env python3
"""Plot candidate-gene microsynteny using GENESPACE combBed and gene_by_species.

The x-axis is gene order, NOT genomic distance. Orthogroup colors show shared
GENESPACE globOG membership, NOT confirmed one-to-one orthology.

Optional GFF3 supplies strand and gene symbols. Without GFF3, genes are drawn
as rectangles to avoid inventing orientation. Requires Python >= 3.7 and matplotlib.

New in v4:
- custom species order via --species-file
- true tree panel from Newick via --tree-file (requires Biopython)
- reference first via topology-preserving node rotation
- species labels can be color-coded by sexual system via --sex-system-file
"""
import argparse
import csv
import gzip
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon, Rectangle, FancyBboxPatch


@dataclass
class Gene:
    species: str
    gene_id: str
    chrom: str
    start: int
    end: int
    glob_og: str
    syn_og: str
    rank: int = -1


def die(message):
    raise SystemExit('ERROR: ' + message)


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')


def read_tsv(path, required):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            die('%s: missing columns %s' % (path, ', '.join(sorted(missing))))
        return list(reader)


def write_tsv(path, rows, columns):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def load_combbed(path):
    rows = read_tsv(path, {'genome', 'id', 'chr', 'start', 'end', 'globOG', 'og'})
    by_chrom = defaultdict(list)
    for r in rows:
        try:
            a, b = int(r['start']), int(r['end'])
        except ValueError:
            continue
        by_chrom[(r['genome'], r['chr'])].append(
            Gene(r['genome'], r['id'], r['chr'], a, b,
                 r['globOG'], r['og']))
    by_id = {}
    for key, genes in by_chrom.items():
        genes.sort(key=lambda g: (g.start, g.end, g.gene_id))
        for i, g in enumerate(genes):
            g.rank = i
            k = g.species, g.gene_id
            if k in by_id:
                die('Duplicate combBed gene ID: %s %s' % k)
            by_id[k] = g
    return by_chrom, by_id


def load_scan(path, gene):
    needed = {'gene', 'species', 'reference_gene_id', 'reference_globOG',
              'reference_synOG', 'local_chr', 'local_gene_ids',
              'family_gene_ids', 'flank_orientation',
              'flank_interval_start_rank', 'flank_interval_end_rank',
              'flank_assessment', 'status'}
    rows = read_tsv(path, needed)
    selected = [r for r in rows if r['gene'] == gene]
    if not selected:
        die('Gene %r not in %s' % (gene, path))
    result = {}
    for row in selected:
        sp = row['species']
        if sp in result:
            die('Duplicate gene/species: %s/%s' % (gene, sp))
        result[sp] = row
    return result


def ordered_unique(items):
    return list(dict.fromkeys(items))


def get_species_and_tree(scan, reference, species_file=None, tree_file=None):
    """Select rows and, if requested, prune/rotate a true Newick topology."""
    if reference not in scan:
        die('Reference %s has no scan row' % reference)
    if species_file:
        with open(species_file) as handle:
            listed = [s.strip() for s in handle if s.strip() and not s.lstrip().startswith('#')]
        if len(listed) != len(set(listed)):
            die('Duplicate names in --species-file')
        missing = [s for s in listed if s not in scan]
        if missing:
            die('Species not in gene_by_species: %s' % ', '.join(missing))
    else:
        listed = list(scan)
    listed = [reference] + [s for s in listed if s != reference]
    if not tree_file:
        return listed, None
    try:
        from Bio import Phylo
    except ImportError:
        die('Tree mode requires Biopython: conda install -c conda-forge biopython')
    try:
        with open_text(tree_file) as handle:
            tree = Phylo.read(handle, 'newick')
    except Exception as exc:
        die('Cannot read a single Newick tree from %s: %s' % (tree_file, exc))
    names = set()
    for tip in tree.get_terminals():
        if not tip.name:
            die('Unnamed terminal in Newick tree')
        # Accept scientific names containing spaces, as well as exact IDs.
        if tip.name not in scan:
            tip.name = tip.name.replace(' ', '_')
        if tip.name in names:
            die('Duplicate terminal name in Newick tree: %s' % tip.name)
        names.add(tip.name)
    wanted = set(listed)
    missing = wanted - names
    if missing:
        die('Missing plotted species in Newick tree: %s. Check names or use --species-file to select a subset.' % ', '.join(sorted(missing)))

    # Keep only the induced subtree. Collapse unary nodes after pruning.
    # Branch lengths are not used in the topology-only drawing.
    def restrict(clade):
        if clade.is_terminal():
            return clade if clade.name in wanted else None
        children = [restrict(child) for child in clade.clades]
        children = [child for child in children if child is not None]
        if not children:
            return None
        if len(children) == 1:
            return children[0]
        clade.clades = children
        return clade

    tree.root = restrict(tree.root)

    # Move the reference-containing child to the first position at each
    # ancestral node. No rerooting, branch swapping, or clade reassignment.
    def rotate(clade):
        if clade.is_terminal():
            return clade.name == reference
        for index, child in enumerate(clade.clades):
            if rotate(child):
                clade.clades.insert(0, clade.clades.pop(index))
                return True
        return False

    rotate(tree.root)
    order = [tip.name for tip in tree.get_terminals()]
    if order[0] != reference or set(order) != wanted:
        die('Internal error: tree tips and plotted species do not match')
    return order, tree


def tree_coordinates(tree, species):
    """Return a rectangular cladogram: aligned tips, no branch-length scale."""
    y_by_name = {sp: len(species)-1-i for i, sp in enumerate(species)}
    depth, y = {}, {}

    def visit(node):
        if node.is_terminal():
            depth[node] = 0
            y[node] = y_by_name[node.name]
        else:
            for child in node.clades:
                visit(child)
            depth[node] = 1 + max(depth[ch] for ch in node.clades)
            ys = [y[ch] for ch in node.clades]
            y[node] = (min(ys) + max(ys)) / 2

    visit(tree.root)
    maximum = depth[tree.root]
    x = {node: maximum-d for node,d in depth.items()}
    return x, y


def draw_tree(ax, tree, species):
    x, y = tree_coordinates(tree, species)
    for node in tree.get_nonterminals():
        child_y = [y[ch] for ch in node.clades]
        ax.plot([x[node], x[node]], [min(child_y), max(child_y)], color='#475569', lw=1.05, solid_capstyle='butt')
        for ch in node.clades:
            ax.plot([x[node], x[ch]], [y[ch], y[ch]], color='#475569', lw=1.05, solid_capstyle='butt')
    maximum = max(x.values())
    ax.set_xlim(-0.2, max(0.2, maximum+0.12))
    ax.axis('off')


def region_for_species(sp, row, reference, by_chrom, by_id, flank):
    """Return a genomic window or None for an unsupported target locus."""
    if sp == reference:
        refid = row['reference_gene_id']
        focal = by_id.get((sp, refid))
        if not focal:
            die('Reference ID %s not found in combBed for %s' % (refid, sp))
        chrom, center, orientation = focal.chrom, float(focal.rank), 1
        family_local = {refid}
        kind = 'reference'
        limits = None
    elif row['flank_assessment'] == 'supported':
        chrom = row['local_chr']
        if (sp, chrom) not in by_chrom:
            die('Supported locus %s/%s has missing chromosome %s' % (row['gene'], sp, chrom))
        try:
            lo = int(row['flank_interval_start_rank'])
            hi = int(row['flank_interval_end_rank'])
        except ValueError:
            die('Bad flank interval ranks for %s/%s' % (row['gene'], sp))
        genes = by_chrom[(sp, chrom)]
        if not (0 <= lo < hi < len(genes)):
            die('Out-of-range interval ranks for %s/%s' % (row['gene'], sp))
        family_local = set(filter(None, row['local_gene_ids'].split(',')))
        invalid = [g for g in family_local
                   if (sp, g) not in by_id or by_id[(sp, g)].chrom != chrom
                   or not lo < by_id[(sp, g)].rank < hi]
        if invalid:
            die('Invalid local_gene_ids for %s/%s: %s' %
                (row['gene'], sp, ','.join(invalid)))
        if family_local:
            center = sum(by_id[(sp, g)].rank for g in family_local) / len(family_local)
            kind = 'local_family'
        else:
            center = (lo + hi) / 2.0
            kind = 'inferred_gap'
        if row['flank_orientation'] not in ('+', '-'):
            die('Missing supported flank orientation for %s/%s' % (row['gene'], sp))
        orientation = 1 if row['flank_orientation'] == '+' else -1
        limits = lo, hi
    elif row['status'] == 'FAMILY_FOUND_LOCUS_UNRESOLVED':
        family_ids = [v.strip() for v in row['family_gene_ids'].split(',') if v.strip()]
        if len(family_ids) != 1:
            return None
        anchor = by_id.get((sp, family_ids[0]))
        if anchor is None or anchor.glob_og != row['reference_globOG']:
            return None
        chrom, center = anchor.chrom, float(anchor.rank)
        orientation = -1 if row['flank_orientation'] == '-' else 1
        family_local = {anchor.gene_id}
        kind = 'family_fallback'
        limits = None
    else:
        return None
    genes = by_chrom[(sp, chrom)]
    window = [(g, orientation * (g.rank - center)) for g in genes
              if abs(g.rank - center) <= flank + 0.51]
    return dict(species=sp, chrom=chrom, orientation=orientation, center=center,
                family_local=family_local, window=window, kind=kind, limits=limits)


def attrs(raw):
    return {key: unquote(val) for part in raw.split(';') if '=' in part
            for key, val in [part.split('=', 1)]}


def gff_file(gff_dir, species):
    candidates = [p for p in gff_dir.iterdir()
                  if p.is_file() and p.name.startswith(species + '.') and
                  p.name.lower().endswith(('.gff', '.gff3', '.gff.gz', '.gff3.gz'))]
    if len(candidates) > 1:
        die('Multiple GFF files for %s: %s' %
            (species, ', '.join(p.name for p in candidates)))
    return candidates[0] if candidates else None


def gff_annotations(path, sp, wanted):
    parents, info, proteins = {}, {}, {}
    wanted_alias = set(wanted)
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            seqid, _, typ, _, _, _, strand, _, raw = fields
            if typ not in ('gene', 'mRNA', 'transcript', 'CDS'):
                continue
            a = attrs(raw)
            fid = a.get('ID', '')
            if typ in ('gene', 'mRNA', 'transcript') and fid:
                parents[fid] = a.get('Parent', '').split(',')[0]
                info[fid] = (strand, a.get('gene') or a.get('gene_name') or a.get('Name') or '', seqid)
            if typ == 'CDS':
                prot = a.get('protein_id', '')
                if not prot and fid.startswith('cds-'):
                    prot = fid[4:]
                if prot in wanted_alias:
                    proteins[prot] = (a.get('Parent', '').split(',')[0], strand, seqid)
    annotated = {}
    for pid in wanted:
        node = pid if pid in info else None
        initial_strand, initial_seqid = '.', ''
        if not node and pid in proteins:
            node, initial_strand, initial_seqid = proteins[pid]
        if not node:
            for prefix in ('gene-', 'rna-'):
                if prefix + pid in info:
                    node = prefix + pid
                    break
        label = ''
        strand, seqid = initial_strand, initial_seqid
        visited = set()
        while node and node not in visited:
            visited.add(node)
            if node in info:
                s, name, chrom = info[node]
                if strand not in ('+', '-') and s in ('+', '-'):
                    strand = s
                if not seqid:
                    seqid = chrom
                if name and not label:
                    label = name
                if node.startswith('gene-') and name:
                    label = name
            node = parents.get(node, '')
        if strand in ('+', '-'):
            annotated[pid] = (strand, label or pid, seqid)
    return annotated


def draw_gene(ax, x, y, color, strand, height=0.32):
    left, right = x - 0.39, x + 0.39
    bottom, top = y - height/2, y + height/2
    if strand not in ('+', '-'):
        patch = FancyBboxPatch((left, bottom), right-left, height,
                               boxstyle='round,pad=0.015,rounding_size=0.035',
                               linewidth=0.45, edgecolor='#475569', facecolor=color,
                               zorder=4)
    else:
        head = 0.16
        if strand == '+':
            vertices = [(left,bottom),(right-head,bottom),(right,y),
                        (right-head,top),(left,top)]
        else:
            vertices = [(right,bottom),(left+head,bottom),(left,y),
                        (left+head,top),(right,top)]
        patch = Polygon(vertices, closed=True, linewidth=0.45,
                        edgecolor='#475569', facecolor=color, zorder=4)
    ax.add_patch(patch)


def elsewhere_text(sp, row, region, by_id):
    if row['status'] != 'FAMILY_MEMBERS_ONLY_ELSEWHERE':
        return ''
    fam = [x for x in row['family_gene_ids'].split(',') if x]
    desc = []
    for gid in fam:
        g = by_id.get((sp, gid))
        if g is not None:
            desc.append('%s:%s' % (g.chrom, g.gene_id))
    if not desc:
        return 'Elsewhere: IDs not found'
    if len(desc) > 2:
        return 'Elsewhere: ' + '; '.join(desc[:2]) + ' (+%d)' % (len(desc)-2)
    return 'Elsewhere: ' + '; '.join(desc)


def load_sex_systems(path):
    if not path:
        return {}
    mapping = {}
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = [x for x in line.replace(',', '\t').split('\t') if x != '']
            if len(parts) < 2:
                die('Bad line in --sex-system-file: %r' % raw.rstrip('\n'))
            species = parts[0].strip()
            system = parts[1].strip()
            mapping[species] = system
    return mapping


def classify_sex_system(label, mode='broad'):
    if not label:
        return ('Unknown', '#64748B')
    s = label.strip().upper().replace(' ', '')
    if mode == 'exact':
        return (label, None)
    if s == 'G':
        return ('G', '#334155')
    if s == 'PG':
        return ('PG', '#BE185D')
    if s == 'PA':
        return ('PA', '#0F766E')
    if '/' in s or '?' in s:
        return ('Mixed/uncertain', '#B45309')
    return (label, '#64748B')


def exact_palette(labels):
    labels = ordered_unique([x for x in labels if x])
    palette = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors)
    return {lab: palette[i % len(palette)] for i, lab in enumerate(labels)}


def run(args):
    cb = Path(args.combbed).expanduser()
    scanpath = Path(args.gene_by_species).expanduser()
    if not cb.is_file() or not scanpath.is_file():
        die('Both --combbed and --gene-by-species must exist')
    scan = load_scan(scanpath, args.gene)
    species, tree = get_species_and_tree(scan, args.reference, args.species_file, args.tree_file)
    by_chrom, by_id = load_combbed(cb)
    sex_systems = load_sex_systems(args.sex_system_file)

    regions = {}
    for sp in species:
        regions[sp] = region_for_species(sp, scan[sp], args.reference, by_chrom, by_id, args.flank)
    refregion = regions[args.reference]
    focal_og = scan[args.reference]['reference_globOG']
    if not focal_og or focal_og in ('NA', '.', 'None'):
        print('WARNING: reference has no valid globOG; only reference focal ID is colored red', file=sys.stderr)
    ref_flanks = []
    for g, x in refregion['window']:
        if (g.gene_id != scan[args.reference]['reference_gene_id'] and g.glob_og and g.glob_og not in ('NA', '.', 'None')):
            if g.glob_og not in ref_flanks and g.glob_og != focal_og:
                ref_flanks.append(g.glob_og)
    palette = list(dict.fromkeys(tuple(c) for name in ('tab20', 'tab20b', 'tab20c') for c in plt.get_cmap(name).colors))
    if len(ref_flanks) > len(palette):
        palette = [plt.get_cmap('hsv')(i / len(ref_flanks)) for i in range(len(ref_flanks))]
    color_of = {og: palette[i] for i, og in enumerate(ref_flanks)}
    fallback = '#D4DBE4'
    focal_color = '#D65749'
    provisional_color = '#D99A3C'

    sex_color_lookup = {}
    if args.sex_system_file:
        if args.sex_color_mode == 'exact':
            sex_color_lookup = exact_palette(sorted(set(sex_systems.values()) | {'Unknown'}))
            sex_color_lookup['Unknown'] = '#64748B'

    annotations = {}
    if args.gff_dir:
        gffdir = Path(args.gff_dir).expanduser()
        if not gffdir.is_dir():
            die('--gff-dir is not a directory: %s' % gffdir)
        for sp in species:
            if not regions[sp]:
                continue
            file = gff_file(gffdir, sp)
            if file is None:
                print('WARNING: No GFF for %s; drawing rectangles' % sp, file=sys.stderr)
                continue
            wanted = {g.gene_id for g, _ in regions[sp]['window']}
            annotations[sp] = gff_annotations(file, sp, wanted)

    n = len(species)
    if tree is not None:
        # Separate axes share identical physical row height and y limits.
        # Widths in inches: tree / names / genes / status information.
        gene_width = max(7.5, (2*args.flank+1)*0.29)
        left_margin, tree_width, name_width, status_width = 0.35, 2.25, 2.9, 4.25
        width = left_margin + tree_width + name_width + gene_width + status_width + 0.3
        # Extra vertical space for rotated reference labels and legends.
        bottom_inches, top_inches = 2.05, 1.1
        row_inches = 0.95 if args.labels == 'all' else 0.58
        track_height = max(2.8, n*row_inches + 0.8)
        height = bottom_inches + track_height + top_inches
        fig = plt.figure(figsize=(width, height))
        bottom, panel_height = bottom_inches/height, track_height/height
        synteny_left = (left_margin+tree_width+name_width)/width
        ax = fig.add_axes([synteny_left, bottom, gene_width/width, panel_height])
        tree_ax = fig.add_axes([left_margin/width, bottom, tree_width/width, panel_height], sharey=ax)
        label_ax = fig.add_axes([(left_margin+tree_width)/width, bottom, name_width/width, panel_height], sharey=ax)
        label_ax.set_xlim(0,1)
        label_ax.axis('off')
        draw_tree(tree_ax, tree, species)
    else:
        height = max(4.6, n * 0.72 + 3.3)
        fig, ax = plt.subplots(figsize=(max(14, args.flank * 1.5 + 8), height))
        fig.subplots_adjust(left=0.30, right=0.75, bottom=0.25, top=0.88)
        label_ax = None
    ax.set_xlim(-args.flank-0.7, args.flank+0.7)
    ax.set_ylim(-0.65, n+1.25)
    ax.set_yticks([])
    for side in ('left', 'right', 'top'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#CBD5E1')
    ax.tick_params(axis='x', colors='#475569', labelsize=8)
    ax.set_xticks(list(range(-args.flank, args.flank + 1, max(1, args.flank//4))))
    ax.set_xlabel('Relative annotated-gene order (not genomic distance)', fontsize=10)
    ax.axvline(0, color='#D65749', linewidth=0.8, linestyle=':', alpha=0.45, zorder=0)
    if tree is not None:
        # Inherit input rooting after pruning; tips are aligned, without a branch-length scale.
        fig.suptitle('Microsynteny of %s | %s reference' % (args.gene, args.reference.replace('_', ' ')), fontsize=14, weight='bold', y=1-0.20/height)
        fig.text((left_margin+tree_width/2)/width, 1-0.78/height, 'Species topology', fontsize=9, weight='bold', ha='center')
        fig.text((left_margin+tree_width+0.12)/width, 1-0.78/height, 'Species / sexual system', fontsize=9, weight='bold')
        fig.text((synteny_left+gene_width/width/2), 1-0.78/height, 'Local gene order', fontsize=9, weight='bold', ha='center')
        fig.text(synteny_left+gene_width/width+0.07/width, 1-0.78/height, 'Locus / assessment', fontsize=9, weight='bold')
    else:
        subtitle = 'species order: custom file' if args.species_file else 'species order: gene_by_species.tsv'
        ax.set_title('Microsynteny of %s (%s reference)\n%s' % (args.gene, args.reference.replace('_', ' '), subtitle), fontsize=14, weight='bold', pad=18)

    if args.links:
        for ix in range(n-1):
            sp1, sp2 = species[ix], species[ix+1]
            r1, r2 = regions[sp1], regions[sp2]
            if not r1 or not r2:
                continue
            by_og1 = defaultdict(list)
            by_og2 = defaultdict(list)
            for g, x in r1['window']:
                by_og1[g.glob_og].append(x)
            for g, x in r2['window']:
                by_og2[g.glob_og].append(x)
            for og in ref_flanks:
                if len(by_og1[og]) == len(by_og2[og]) == 1:
                    provisional = ('family_fallback' in (r1['kind'], r2['kind']))
                    ax.plot([by_og1[og][0], by_og2[og][0]], [n-1-ix, n-2-ix], color=color_of[og],
                            linestyle='--' if provisional else '-', linewidth=1.15,
                            alpha=0.18 if provisional else 0.3, zorder=1)

    gene_rows, summary_rows = [], []
    broad_seen = set()
    for ix, sp in enumerate(species):
        y = n-1-ix
        row, region = scan[sp], regions[sp]
        sex_raw = sex_systems.get(sp, '')
        if args.sex_system_file:
            if args.sex_color_mode == 'exact':
                sex_group = sex_raw if sex_raw else 'Unknown'
                sex_color = sex_color_lookup.get(sex_group, '#64748B')
            else:
                sex_group, sex_color = classify_sex_system(sex_raw, 'broad')
                broad_seen.add((sex_group, sex_color))
        else:
            sex_group, sex_color = '', '#111827'
        label = sp.replace('_', ' ')
        if args.show_sex_label and sex_raw:
            label = '%s [%s]' % (label, sex_raw)
        if label_ax is not None:
            label_ax.text(0.045, y, label, ha='left', va='center', fontsize=9, fontstyle='italic',
                          fontweight='bold' if sp == args.reference else 'normal', color=sex_color, clip_on=False)
        else:
            ax.text(-0.02, y, label, ha='right', va='center', fontsize=9, fontstyle='italic',
                    color=sex_color, transform=ax.get_yaxis_transform(), clip_on=False)
        if region is None:
            ax.text(0, y, 'Unresolved locus: no reliable two-sided flank support', fontsize=8.5, color='#7A8798', ha='center', va='center')
            warning = elsewhere_text(sp,row,region,by_id)
            ax.text(1.015, y+0.12, 'UNRESOLVED', fontsize=8, color='#A36B19', ha='left', va='center', transform=ax.get_yaxis_transform())
            if warning:
                ax.text(1.015, y-0.15, warning, fontsize=6.6, color='#B05B30', ha='left', va='center', transform=ax.get_yaxis_transform())
            summary_rows.append(dict(gene=args.gene, species=sp, status=row['status'], displayed='no', chromosome='', anchor_kind='unresolved', orientation='', family_elsewhere=elsewhere_text(sp,row,region,by_id), sex_system=sex_raw, sex_group=sex_group))
            continue
        ax.hlines(y, -args.flank-0.46, args.flank+0.46, linewidth=0.85, color='#AFBCCB', zorder=0)
        if region['kind'] == 'inferred_gap':
            lo, hi = region['limits']
            left = region['orientation'] * (lo-region['center'])
            right = region['orientation'] * (hi-region['center'])
            a, b = sorted((left,right))
            a, b = max(a,-args.flank-0.55), min(b,args.flank+0.55)
            if b > a:
                ax.add_patch(Rectangle((a,y-0.21),b-a,0.42, fill=False, linestyle='--', linewidth=1.1, edgecolor=focal_color, alpha=0.8, zorder=2))
        for g,x in region['window']:
            is_focal = (sp == args.reference and g.gene_id == row['reference_gene_id']) or (g.gene_id in region['family_local'] and sp != args.reference)
            color = provisional_color if (is_focal and region['kind'] == 'family_fallback') else (focal_color if is_focal else color_of.get(g.glob_og, fallback))
            raw_strand, label2, source_chrom = annotations.get(sp,{}).get(g.gene_id, ('.',g.gene_id,''))
            if source_chrom and source_chrom != g.chrom:
                raw_strand, label2 = '.', g.gene_id
            shown_strand = raw_strand
            if region['orientation'] == -1:
                shown_strand = {'+':'-','-':'+'}.get(raw_strand,'.')
            draw_gene(ax,x,y,color,shown_strand)
            gene_rows.append(dict(gene=args.gene,species=sp,chromosome=g.chrom,gene_id=g.gene_id,
                                  combbed_start=g.start,combbed_end=g.end,rank=g.rank,
                                  plotted_x=(0.0 if abs(x) < 1e-9 else round(x,3)),globOG=g.glob_og,
                                  synOG=g.syn_og,original_strand=raw_strand,displayed_strand=shown_strand,
                                  gene_label=label2,focal_family_in_expected_locus=int(is_focal and region['kind'] != 'family_fallback'),
                                  status=row['status'],sex_system=sex_raw,sex_group=sex_group))
            if ((sp==args.reference and args.labels in ('reference','all')) or (sp!=args.reference and args.labels == 'all')):
                if is_focal or g.glob_og in color_of:
                    name = args.gene if is_focal else label2
                    ax.text(x-0.25,y+0.24,name,fontsize=5.8,rotation=55,ha='left',va='bottom',color='#334155',clip_on=False)
        if region['kind'] == 'family_fallback':
            at = '%s | FAMILY HIT (synteny unresolved)' % region['chrom']
        else:
            at = '%s | %s' % (region['chrom'], row['status'].replace('FAMILY_MEMBERS_ONLY_ELSEWHERE','ELSEWHERE').replace('SYNTENIC_FAMILY_MEMBER','SYNTENIC').replace('NO_ANNOTATED_MEMBER_IN_CONSERVED_LOCUS','NO_ANNOTATED_MEMBER'))
        ax.text(1.015, y+0.12, at, fontsize=7.4, color='#334155', ha='left', va='center', transform=ax.get_yaxis_transform())
        off = elsewhere_text(sp, row, region, by_id)
        if off:
            ax.text(1.015,y-0.15,off,fontsize=6.6,color='#B05B30',ha='left',va='center',transform=ax.get_yaxis_transform())
        summary_rows.append(dict(gene=args.gene, species=sp, status=row['status'], displayed='yes', chromosome=region['chrom'], anchor_kind=region['kind'], orientation=('+' if region['orientation']==1 else '-'), family_elsewhere=off, sex_system=sex_raw, sex_group=sex_group))

    legend1 = [Patch(facecolor=focal_color,edgecolor='#475569',label='Candidate family at expected locus'),
               Patch(facecolor=provisional_color,edgecolor='#475569',label='Unique family hit (synteny unresolved)'),
               Patch(facecolor=palette[0],edgecolor='#475569',label='Reference-flank globOG (same color = same OG)'),
               Patch(facecolor=fallback,edgecolor='#475569',label='Other annotated genes'),
               Patch(facecolor='none',edgecolor=focal_color,linestyle='--',label='Flank-inferred interval; no local family gene')]
    fig.legend(handles=legend1,loc='lower center',bbox_to_anchor=(0.5,0.57/height if tree is not None else 0.045), ncol=3,fontsize=8,frameon=False)

    if args.sex_system_file:
        if args.sex_color_mode == 'exact':
            handles = []
            for k in ordered_unique([sex_systems.get(sp, 'Unknown') for sp in species]):
                handles.append(Patch(facecolor=sex_color_lookup.get(k, '#64748B'), edgecolor='none', label=k))
        else:
            ordered_groups = [('G', '#334155'), ('PG', '#BE185D'), ('PA', '#0F766E'), ('Mixed/uncertain', '#B45309'), ('Unknown', '#64748B')]
            handles = [Patch(facecolor=c, edgecolor='none', label=g) for g,c in ordered_groups if any(g == seen[0] for seen in broad_seen)]
        if tree is not None:
            fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5,0.12/height), ncol=min(8,max(1,len(handles))), title='Sexual system', fontsize=8, title_fontsize=8, frameon=False)
        else:
            fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.012,0.90), title='Sexual system', fontsize=8, title_fontsize=8, frameon=False)

    if tree is not None:
        fig.text(0.5, 1.21/height, 'Tree shows induced topology (input rooting; rotated nodes; no branch-length scale). Gene arrows require matching GFF3 strand data.', ha='center', fontsize=8, color='#64748B')
        fig.text(0.5, 1.03/height, 'Orange family hits and dashed links are provisional; a reversed display is not proof of a genomic inversion.', ha='center', fontsize=8, color='#9A691B')
    else:
        fig.text(0.30,0.165,'Strands shown as arrows only when supported by matching GFF3 IDs; reversed loci are displayed in reference orientation.', fontsize=8,color='#64748B')
        fig.text(0.30,0.142,'Orange family-anchored rows are provisional; orientation from failed flank assessment is not validated.', fontsize=8,color='#9A691B')
        if args.sex_system_file:
            fig.text(0.30,0.119,'Species labels are color-coded by sexual system (%s mode).' % args.sex_color_mode, fontsize=8, color='#64748B')
    base = Path(args.output_prefix).expanduser()
    base.parent.mkdir(parents=True,exist_ok=True)
    for ext in ('svg','pdf','png'):
        fig.savefig(str(base)+'.'+ext,dpi=250,facecolor='white')
    plt.close(fig)
    if tree is not None:
        from Bio import Phylo
        # Topology-only Newick, matching the displayed order exactly.
        Phylo.write(tree, str(base)+'.display_tree.nwk', 'newick', plain=True)
    with open(str(base)+'.species_order.txt','w') as handle:
        handle.write('\n'.join(species)+'\n')
    fields1 = ['gene','species','chromosome','gene_id','combbed_start','combbed_end','rank','plotted_x','globOG','synOG','original_strand','displayed_strand','gene_label','focal_family_in_expected_locus','status','sex_system','sex_group']
    fields2 = ['gene','species','status','displayed','chromosome','anchor_kind','orientation','family_elsewhere','sex_system','sex_group']
    write_tsv(str(base)+'.genes.tsv',gene_rows,fields1)
    write_tsv(str(base)+'.regions.tsv',summary_rows,fields2)
    print('Saved: '+', '.join(str(base)+'.'+x for x in ('svg','pdf','png','genes.tsv','regions.tsv')))
    if not args.gff_dir:
        print('NOTE: GFF3 not supplied. Rectangles are used (strand is unknown).')
    print('NOTE: globOG family matching and flank intervals are candidates, not proven orthology/loss.')
    print('NOTE: FAMILY HIT (synteny unresolved) rows are centered on a unique family member, not a validated locus.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--combbed', required=True,help='GENESPACE results/combBed.txt')
    p.add_argument('--gene-by-species',required=True,help='Output from sparid_sexgene_scan.py')
    p.add_argument('--gene',required=True,help='Candidate gene symbol, e.g. dmrt1')
    p.add_argument('--reference',default='Pagrus_major')
    p.add_argument('--species-file',help='One species per line: order without tree; subset selection with tree')
    p.add_argument('--tree-file',help='Newick file: draw rooted topology at left (requires Biopython)')
    p.add_argument('--gff-dir',help='Optional directory with Species_name*.gff[3][.gz] files')
    p.add_argument('--sex-system-file',help='Two-column TSV/CSV: species<TAB>sexual_system')
    p.add_argument('--sex-color-mode',choices=('broad','exact'),default='broad', help='Broad: G/PG/PA/Mixed-uncertain. Exact: each raw label gets its own color')
    p.add_argument('--show-sex-label',action='store_true',help='Append raw sex-system code after species names')
    p.add_argument('--flank',type=int,default=10,help='Genes per side to visualize')
    p.add_argument('--labels',choices=('none','reference','all'),default='reference')
    p.add_argument('--links',action='store_true',help='Connect unique visible flank globOG members in adjacent rows')
    p.add_argument('--output-prefix',required=True,help='Output path without extension')
    args = p.parse_args()
    if args.flank < 2 or args.flank > 50:
        p.error('--flank must be between 2 and 50')
    run(args)


if __name__ == '__main__':
    main()
