#!/usr/bin/env perl
use strict;
use warnings;

# ============================================================
# go_enrichment_assay.pl  (ALL-IN-ONE)
# ============================================================
# What this single script contains:
#   1) A lightweight DBM helper ("DBMUtil") that replaces DBM.pm
#      - loads MLDBM / Storable(.nstore) hashes
#      - supports a simple KEY->PATH mapping file (TSV) instead of hard-coded paths
#
#   2) GO enrichment engine ("GOPhyper") adapted from GO_phyper_coffee.pm
#      - uses Statistics::R::phyper to compute hypergeometric p-values
#      - performs FDR (Benjamini-Hochberg)
#      - reads target gene list files and writes enrichment table
#
#   3) A generalized main runner that replaces 3_each_Node_GO_phyper.pl
#      - reads a "targets TSV" describing each node + which species DBs to combine
#      - resolves all DBM paths via your mapping TSV (portable for GitHub)
#
# ------------------------------------------------------------
# REQUIRED INPUT FORMAT (portable, no hard-coded windows paths)
#
# (A) Mapping TSV (KEY -> PATH), e.g. mldbm_paths.tsv
#     key<TAB>path
#     arabica.go.id<TAB>/path/to/ArabicaMrna2GoidDBM
#     arabica.GO.id.Ancestral<TAB>/path/to/ArabicaMrna2AncestralGoidDBM
#     arabica.GO.Count.ancestral<TAB>/path/to/ArabicaAncestralGOCountDBM
#     eug.go.id<TAB>/path/to/EugenoidMrna2GoidDBM
#     ...
#     GOID2TermDBM<TAB>/path/to/GOID2TermDBM
#
# (B) Targets TSV, e.g. targets.tsv
#     node_id<TAB>mldbm_prefixes<TAB>gene_file_regex<TAB>out_csv
#     27_expansion<TAB>arabica.eug.robusta<TAB>^27_expansion_geneList_\.tab$<TAB>27_expansion_GO_NoAncestral.csv
#
# (C) Gene list files (in --gene-dir), matching gene_file_regex
#     The gene list file is assumed to contain gene IDs in the FIRST column (TSV).
#     (extra columns are ignored)
#
# ------------------------------------------------------------
# USAGE EXAMPLE
#
#   perl go_enrichment_single.pl \
#     --map mldbm_paths.tsv \
#     --targets targets.tsv \
#     --gene-dir . \
#     --pvalue 0.01 \
#     --aspect P \
#     --bg genome \
#     --use-no-ancestral-target-go 1
#
# If you want ancestral GO for TARGET genes too:
#     --use-no-ancestral-target-go 0
#
# Notes:
# - Background GO count should usually be the "ancestral GO count" DB(s),
#   even if target GO uses non-ancestral (this matches your notes).
# - This script expects you have R + Statistics::R working.
# ============================================================

use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);

# ----------------------------
# CLI
# ----------------------------
my $map_tsv     = '';
my $targets_tsv = '';
my $gene_dir    = '.';
my $pvalue_cut  = 0.01;
my $fdr_cut     = 0.05;
my $aspect      = 'P';          # P/F/C
my $bg          = 'genome';     # genome OR path to bg gene list file
my $use_no_anc_target_go = 1;   # 1 = target GO uses .go.id ; 0 = target GO uses .GO.id.Ancestral
my $verbose     = 1;

my $help = 0;

GetOptions(
  'map=s'      => \$map_tsv,
  'targets=s'  => \$targets_tsv,
  'gene-dir=s' => \$gene_dir,

  'pvalue=f'   => \$pvalue_cut,
  'fdr=f'      => \$fdr_cut,
  'aspect=s'   => \$aspect,
  'bg=s'       => \$bg,

  'use-no-ancestral-target-go=i' => \$use_no_anc_target_go,
  'verbose=i'  => \$verbose,

  'help!'      => \$help,
) or die usage();

if ($help) { print usage(); exit 0; }
die "[ERROR] --map is required\n" unless $map_tsv;
die "[ERROR] --targets is required\n" unless $targets_tsv;

$gene_dir = abs_path($gene_dir);
$map_tsv = abs_path($map_tsv);
$targets_tsv = abs_path($targets_tsv);

die "[ERROR] gene-dir not found: $gene_dir\n" unless -d $gene_dir;
die "[ERROR] map TSV not found: $map_tsv\n" unless -f $map_tsv;
die "[ERROR] targets TSV not found: $targets_tsv\n" unless -f $targets_tsv;

# Load KEY->PATH mapping
my $MAP = DBMUtil::read_map_tsv($map_tsv);

# GOID2TermDBM must be provided in map (key "GOID2TermDBM")
die "[ERROR] Mapping TSV missing key: GOID2TermDBM\n" unless exists $MAP->{GOID2TermDBM};

# Read targets and run
my $targets = read_targets_tsv($targets_tsv);

for my $t (@$targets) {
  my ($node_id, $prefixes_dot, $gene_regex, $out_csv) = @$t;

  print "\n=== NODE: $node_id ===\n" if $verbose;

  # prefixes like arabica.eug.robusta
  my @prefixes = split(/\./, $prefixes_dot);

  # Resolve DBM keys for this node:
  # - bg_go_count:  <prefix>.GO.Count.ancestral
  # - genome_total_GO_Ancestral: <prefix>.GO.id.Ancestral
  # - genome_total_GO (No-ancestral target GO): <prefix>.go.id   (note lower-case in your old script)
  #
  my @bg_go_count_dbm = map {
    my $k = "$_.GO.Count.ancestral";
    die "[ERROR] map missing key: $k\n" unless exists $MAP->{$k};
    $MAP->{$k}
  } @prefixes;

  my @genome_anc_dbm = map {
    my $k = "$_.GO.id.Ancestral";
    die "[ERROR] map missing key: $k\n" unless exists $MAP->{$k};
    $MAP->{$k}
  } @prefixes;

  my @genome_noanc_dbm;
  if ($use_no_anc_target_go) {
    @genome_noanc_dbm = map {
      my $k = "$_.go.id";
      die "[ERROR] map missing key: $k\n" unless exists $MAP->{$k};
      $MAP->{$k}
    } @prefixes;
  } else {
    # If target GO uses ancestral too, just reuse the ancestral DBM(s)
    @genome_noanc_dbm = @genome_anc_dbm;
  }

  # Find gene-list files matching this node
  opendir(my $DIR, $gene_dir) or die "[ERROR] cannot open gene-dir: $gene_dir ($!)\n";
  my @files = grep { /$gene_regex/ } readdir($DIR);
  closedir($DIR);

  if (!@files) {
    warn "[WARN] No gene list file matched regex [$gene_regex] in $gene_dir; skipping $node_id\n";
    next;
  }
  print "[INFO] matched files: ", join(", ", @files), "\n" if $verbose;

  my $go = GOPhyper->new(
    dir                     => $gene_dir,
    pvalue                  => $pvalue_cut,
    fdr                     => $fdr_cut,
    aspect                  => $aspect,
    bg                      => $bg,                 # "genome" or file path to bg genes
    bg_go_count             => \@bg_go_count_dbm,    # combined
    goid2term               => $MAP->{GOID2TermDBM},
    genome_total_GO_Ancestral => \@genome_anc_dbm,   # combined
    genome_total_GO         => \@genome_noanc_dbm,   # combined (target GO set)
    verbose                 => $verbose,
  );

  # Run enrichment for matched gene list files
  $go->print_enrichment_score(\@files, $out_csv);
}

print "\n[DONE]\n";
exit 0;

# ============================================================
# TARGETS TSV PARSER
# ============================================================
sub read_targets_tsv {
  my ($path) = @_;
  open(my $IN, "<", $path) or die "[ERROR] cannot open targets TSV: $path ($!)\n";
  my @rows;
  while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*#/;
    my @c = split(/\t/, $line, -1);
    die "[ERROR] targets.tsv must have 4 columns: node_id, mldbm_prefixes, gene_file_regex, out_csv\nLine: $line\n"
      unless @c >= 4;
    push @rows, [ @c[0..3] ];
  }
  close $IN;
  return \@rows;
}

sub usage {
  return <<'USAGE';
Usage:
  perl go_enrichment_assay.pl --map mldbm_paths.tsv --targets targets.tsv [options]

Required:
  --map      TSV with columns: key<TAB>path
  --targets  TSV with columns: node_id<TAB>mldbm_prefixes<TAB>gene_file_regex<TAB>out_csv

Options:
  --gene-dir DIR     Directory containing gene list files (default: .)
  --pvalue FLOAT     P-value cutoff (default: 0.01)
  --fdr FLOAT        FDR cutoff for BH adjustment (default: 0.05)
  --aspect STR       GO aspect: P|F|C (default: P)
  --bg STR           Background: "genome" OR a gene list file path (default: genome)
  --use-no-ancestral-target-go 0|1
                     1 = use "<prefix>.go.id" for target genes (default)
                     0 = use "<prefix>.GO.id.Ancestral" for target genes too
  --verbose 0|1      Print progress (default: 1)

Mapping keys expected (per species prefix):
  <prefix>.go.id
  <prefix>.GO.id.Ancestral
  <prefix>.GO.Count.ancestral
And also:
  GOID2TermDBM

Example mapping TSV:
  arabica.go.id     /path/ArabicaMrna2GoidDBM
  arabica.GO.id.Ancestral /path/ArabicaMrna2AncestralGoidDBM
  arabica.GO.Count.ancestral /path/ArabicaAncestralGOCountDBM
  GOID2TermDBM      /path/GOID2TermDBM

Example targets TSV:
  27_expansion  arabica.eug.robusta  ^27_expansion_geneList_\.tab$  27_expansion_GO_NoAncestral.csv
USAGE
}

# ============================================================
# PACKAGE: DBMUtil   (replacement for DBM.pm in one file)
# ============================================================
{
  package DBMUtil;
  use strict;
  use warnings;
  use Carp qw(confess);
  use Storable qw(retrieve);
  use MLDBM qw(DB_File Storable portable);
  use Fcntl;

  # Read a portable mapping TSV: key<TAB>path
  sub read_map_tsv {
    my ($path) = @_;
    open(my $IN, "<", $path) or die "[ERROR] cannot open map TSV: $path ($!)\n";
    my %m;
    while (my $line = <$IN>) {
      chomp $line;
      next if $line =~ /^\s*$/;
      next if $line =~ /^\s*#/;
      my ($k, $v) = split(/\t/, $line, 2);
      next unless defined $k and defined $v;
      $k =~ s/^\s+|\s+$//g;
      $v =~ s/^\s+|\s+$//g;
      $m{$k} = $v if length($k) and length($v);
    }
    close $IN;
    return \%m;
  }

  # Load either:
  #  - MLDBM directory/file (tie)
  #  - Storable file (.nstore / .store / .storable)
  #
  # Returns a HASHREF (for MLDBM, it returns the tied hashref).
  sub retrieveDbm {
    my ($path_or_hashref) = @_;
    # Allow caller to pass a hashref directly (combined output)
    if (ref($path_or_hashref) eq 'HASH') {
      return $path_or_hashref;
    }
    my $path = $path_or_hashref;

    confess "retrieveDbm: no path given" unless defined $path && length($path);

    if (-f $path && $path =~ /\.(nstore|store|storable)$/i) {
      my $h = retrieve($path);
      confess "retrieveDbm: Storable did not return a HASH ref for $path" unless ref($h) eq 'HASH';
      return $h;
    }

    # MLDBM: often is a file prefix or directory; tie works with DB_File backend
    my %h;
    tie %h, 'MLDBM', $path, O_RDONLY, 0666 or confess "retrieveDbm: cannot tie MLDBM ($path): $!";
    return \%h;
  }

  # Combine multiple DBMs by summing counts (for GO-count DBMs) or merging hashes.
  # For count tables: key => number. For other tables: key => anything (keeps first).
  sub combine_dbm_sum {
    my ($paths_aref) = @_;
    confess "combine_dbm_sum expects ARRAY ref" unless ref($paths_aref) eq 'ARRAY';
    my %combined;
    for my $p (@$paths_aref) {
      my $h = retrieveDbm($p);
      for my $k (keys %$h) {
        if (!exists $combined{$k}) {
          $combined{$k} = $h->{$k};
        } else {
          # sum numeric if possible
          if (defined $combined{$k} && defined $h->{$k} && $combined{$k} =~ /^\d+$/ && $h->{$k} =~ /^\d+$/) {
            $combined{$k} += $h->{$k};
          }
        }
      }
    }
    return \%combined;
  }

  # Combine multiple genome GO annotation DBMs:
  # gene => {P=>[GO..],F=>[],C=>[]}
  # If same gene appears in multiple DBs, later one overwrites (rare for different species).
  sub combine_dbm_merge {
    my ($paths_aref) = @_;
    confess "combine_dbm_merge expects ARRAY ref" unless ref($paths_aref) eq 'ARRAY';
    my %combined;
    for my $p (@$paths_aref) {
      my $h = retrieveDbm($p);
      for my $k (keys %$h) {
        $combined{$k} = $h->{$k};
      }
    }
    return \%combined;
  }
}

# ============================================================
# PACKAGE: GOPhyper  (embedded + generalized from GO_phyper_coffee.pm)
# ============================================================
{
  package GOPhyper;
  use strict;
  use warnings;
  use Carp qw(confess);
  use Cwd qw(abs_path);
  use File::Basename qw(basename);
  use Statistics::R;

  sub new {
    my ($class, %args) = @_;

    my $self = bless {}, $class;

    $self->{_dir}    = $args{dir}    // '.';
    $self->{_pval}   = $args{pvalue} // 0.01;
    $self->{_fdr}    = $args{fdr}    // 0.05;
    $self->{_aspect} = $args{aspect} // 'P';
    $self->{_bg}     = $args{bg}     // 'genome';
    $self->{_verbose}= $args{verbose}// 1;

    confess "goid2term is required" unless $args{goid2term};
    confess "bg_go_count is required" unless $args{bg_go_count};
    confess "genome_total_GO_Ancestral is required" unless $args{genome_total_GO_Ancestral};
    confess "genome_total_GO is required" unless $args{genome_total_GO};

    # bg_go_count can be ARRAY(ref) of MLDBM paths (sum counts)
    my $bg_go_count_tbl =
      (ref($args{bg_go_count}) eq 'ARRAY')
        ? DBMUtil::combine_dbm_sum($args{bg_go_count})
        : DBMUtil::retrieveDbm($args{bg_go_count});

    # goid2term can be single DBM path
    my $goid2term_tbl = DBMUtil::retrieveDbm($args{goid2term});

    # genome_total_GO_Ancestral can be multiple species (merge)
    my $genome_anc_tbl =
      (ref($args{genome_total_GO_Ancestral}) eq 'ARRAY')
        ? DBMUtil::combine_dbm_merge($args{genome_total_GO_Ancestral})
        : DBMUtil::retrieveDbm($args{genome_total_GO_Ancestral});

    # genome_total_GO (for target GO extraction; often no-ancestral)
    my $genome_noanc_tbl =
      (ref($args{genome_total_GO}) eq 'ARRAY')
        ? DBMUtil::combine_dbm_merge($args{genome_total_GO})
        : DBMUtil::retrieveDbm($args{genome_total_GO});

    # Determine background set
    my ($bg_GO_hash, $bg_gene_number);
    if ($self->{_bg} =~ /^genome$/i) {
      $bg_GO_hash = $genome_anc_tbl;
      $bg_gene_number = scalar(keys %$genome_anc_tbl);
    } elsif (-f $self->{_bg}) {
      my $bg_gene_list = $self->make_gene_list_from_file($self->{_bg});
      $bg_GO_hash = $self->slice_hash($genome_anc_tbl, $bg_gene_list);
      $bg_gene_number = scalar(@$bg_gene_list);
    } else {
      confess "bg must be 'genome' or an existing file path; got: $self->{_bg}";
    }

    $self->{_bg_go_count}   = $bg_go_count_tbl;
    $self->{_goid2term}     = $goid2term_tbl;
    $self->{_genome_anc}    = $genome_anc_tbl;
    $self->{_genome_noanc}  = $genome_noanc_tbl;
    $self->{_bg_anno}       = $bg_GO_hash;
    $self->{_bg_total_n}    = $bg_gene_number;

    $self->{_R} = $self->_Rstart($self->{_dir});
    return $self;
  }

  sub _Rstart {
    my ($self, $dir) = @_;
    my $R = Statistics::R->new();
    $dir =~ s/\\/\//g;
    $R->run(qq'setwd("$dir")');
    if ($self->{_verbose}) {
      my $cwd = $R->run('print(getwd())');
      print "\n== R_getwd(): $cwd\n";
    }
    return $R;
  }

  # ----------------------------------------------------------
  # Main: run GO enrichment for each gene list file
  # ----------------------------------------------------------
  sub print_enrichment_score {
    my ($self, $files_aref, $out_csv) = @_;
    confess "files must be ARRAY ref" unless ref($files_aref) eq 'ARRAY';
    confess "out_csv required" unless $out_csv;

    my $R = $self->{_R};

    open(my $OUT, ">", $out_csv) or die "[ERROR] cannot write $out_csv ($!)\n";
    print $OUT join(",", qw(
      NodeFile Aspect GOID Term
      bg_total_genes target_total_genes
      bg_GO_count target_GO_count
      pvalue fdr
    )), "\n";

    for my $file (@$files_aref) {
      my $path = $self->{_dir} . "/" . $file;

      my $target_genes = $self->make_gene_list_from_file($path);
      my %target_set = map { $_ => 1 } @$target_genes;

      my $target_total = scalar(@$target_genes);
      my $bg_total     = $self->{_bg_total_n};

      print "[INFO] $file : bg_total=$bg_total target_total=$target_total\n" if $self->{_verbose};

      # Count GO occurrences in target genes using genome_total_GO (often NO-ancestral)
      my $target_go_counts = $self->count_target_go(\%target_set);

      # For each GOID in target, compute hypergeometric over-representation p-value
      my @rows;
      my @pvals;

      for my $goid (sort keys %$target_go_counts) {
        my $target_k = $target_go_counts->{$goid};             # successes in sample
        my $bg_M     = $self->{_bg_go_count}{$goid} // 0;      # successes in population (background)
        next if $bg_M <= 0;

        # Hypergeometric parameters:
        #   population size = bg_total
        #   number of successes in population = bg_M
        #   sample size = target_total
        #   successes in sample = target_k
        #
        # Over-representation p-value:
        #   P(X >= target_k) = phyper(target_k-1, bg_M, bg_total-bg_M, target_total, lower.tail=FALSE)
        #
        my $p = $self->phyper_righttail($target_k, $bg_M, $bg_total, $target_total);

        push @rows, [$file, $self->{_aspect}, $goid, ($self->{_goid2term}{$goid} // 'NA'),
                     $bg_total, $target_total, $bg_M, $target_k, $p, 'NA'];
        push @pvals, $p;
      }

      # Apply BH FDR
      my $fdrs = bh_adjust(\@pvals);
      for (my $i=0; $i<@rows; $i++) {
        $rows[$i]->[9] = $fdrs->[$i]; # fdr
      }

      # Write filtered results by pvalue cutoff (and optionally fdr cutoff)
      for my $r (@rows) {
        my ($nodefile, $asp, $goid, $term, $bgN, $tN, $bgM, $tk, $p, $fdr) = @$r;
        next if $p > $self->{_pval};
        # If you want to also filter by FDR, uncomment:
        # next if $fdr > $self->{_fdr};

        # Clean commas from term for CSV safety
        $term =~ s/,/ /g;

        print $OUT join(",", map { csv_escape($_) } @$r), "\n";
      }
    }

    close $OUT;
    print "[OK] wrote $out_csv\n" if $self->{_verbose};
  }

  # Count target GO using genome_total_GO table:
  # genome_total_GO: gene => {P=>[GO..],F=>[],C=>[]}
  sub count_target_go {
    my ($self, $target_set) = @_;
    my %cnt;
    my $tbl = $self->{_genome_noanc};

    for my $gene (keys %$target_set) {
      next unless exists $tbl->{$gene};
      my $asp = $self->{_aspect};
      my $goids = $tbl->{$gene}{$asp};
      next unless $goids && ref($goids) eq 'ARRAY';
      for my $go (@$goids) {
        $cnt{$go}++;
      }
    }
    return \%cnt;
  }

  # Right-tail hypergeometric using R phyper
  sub phyper_righttail {
    my ($self, $k, $M, $N, $n) = @_;
    my $R = $self->{_R};

    # phyper(q, m, n, k, lower.tail=FALSE)
    # here: q = k-1, m = M, n = N-M, k = n(sample size)
    my $q = $k - 1;
    $q = 0 if $q < 0;

    # Use R to compute
    $R->set('q', $q);
    $R->set('m', $M);
    $R->set('n', $N - $M);
    $R->set('k', $n);
    my $p = $R->run('as.numeric(phyper(q, m, n, k, lower.tail=FALSE))');
    $p =~ s/^\s+|\s+$//g;
    $p = 1.0 if $p eq '' || $p !~ /^[0-9eE\.\-\+]+$/;
    return $p + 0; # numeric
  }

  # Gene list reader:
  #   - expects gene IDs in first column
  #   - skips header lines starting with '#' or empty lines
  sub make_gene_list_from_file {
    my ($self, $file) = @_;
    open(my $IN, "<", $file) or confess "Cannot open gene list file: $file ($!)";
    my @genes;
    while (my $line = <$IN>) {
      chomp $line;
      next if $line =~ /^\s*$/;
      next if $line =~ /^\s*#/;
      my ($gid) = split(/\t/, $line, 2);
      next unless defined $gid && length $gid;
      $gid =~ s/^\s+|\s+$//g;
      next unless length $gid;
      push @genes, $gid;
    }
    close $IN;
    return \@genes;
  }

  sub slice_hash {
    my ($self, $hashref, $keys_aref) = @_;
    my %s;
    for my $k (@$keys_aref) {
      $s{$k} = $hashref->{$k} if exists $hashref->{$k};
    }
    return \%s;
  }

  # BH adjustment (Benjamini-Hochberg)
  sub bh_adjust {
    my ($pvals) = @_;
    my $n = scalar(@$pvals);
    return [] if $n == 0;

    # keep original indices
    my @idx = sort { $pvals->[$a] <=> $pvals->[$b] } 0..$#$pvals;
    my @q = (0) x $n;

    my $prev = 1;
    for (my $rank=$n; $rank>=1; $rank--) {
      my $i = $idx[$rank-1];
      my $p = $pvals->[$i];
      my $adj = ($p * $n) / $rank;
      $adj = 1 if $adj > 1;
      # enforce monotonicity from largest rank downwards
      $prev = $adj if $adj < $prev;
      $q[$i] = $prev;
    }
    return \@q;
  }

  sub csv_escape {
    my ($v) = @_;
    $v = '' unless defined $v;
    # If it contains comma/quote/newline, wrap in quotes and escape quotes
    if ($v =~ /[,"\n]/) {
      $v =~ s/"/""/g;
      return qq("$v");
    }
    return $v;
  }
}

# end
