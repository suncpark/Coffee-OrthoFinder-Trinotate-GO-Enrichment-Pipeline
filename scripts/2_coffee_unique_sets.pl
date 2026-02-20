#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use List::Util qw(sum);

# ============================================================
# coffee_unique_sets.pl
#   Combine two utilities into ONE script:
#
#   (A) Coffee-unique orthogroups (from familywise count table)
#       Input  : 0_each_ortho_familywiseCount.tsv
#       Output : ara_unique_ortho.tsv, eug_unique_ortho.tsv, rob_unique_ortho.tsv,
#                ara.eug_unique_ortho.tsv, ara.rob_unique_ortho.tsv, eug.rob_unique_ortho.tsv,
#                ara.eug.rob_unique_ortho.tsv
#
#   (B) Coffee-unique unassigned genes (from Orthogroups_UnassignedGenes.tsv)
#       Input  : Orthogroups_UnassignedGenes.tsv
#       Output : arabica_unique_unassigned.tsv, eug_unique_unassigned.tsv, robusta_unique_unassigned.tsv
#
# Notes:
#   - Mode A reproduces your 1_retrieve_Coffee_Unique_Orthogroups.pl logic.
#   - Mode B implements the intended behavior of 2_retrieve_Coffee_Unique_Unassigned.pl
#     (the uploaded file contained notes only).
# ============================================================

my $familywise = '../1_orthologs_familyWise/0_each_ortho_familywiseCount.tsv';
my $unassigned = '../Orthogroups_UnassignedGenes.tsv';
my $outdir     = '.';

# For familywise-count table: default columns are assumed to be:
#   Orthogroup  laminales  solnales  olden  arabica  eug  robusta
# (6 category columns after Orthogroup)
# If your file differs, use --family-cols to define the 6 column names (comma-separated).
my $family_cols = 'laminales,solnales,olden,arabica,eug,robusta';

# For unassigned table, we need to locate coffee columns by header matching.
# You can give exact header names or regex fragments.
my $arabica_col_pat = 'CoffeaArabica|arabica';
my $eug_col_pat     = 'CoffeaEugenioides|eugenioides|eug';
my $rob_col_pat     = 'CoffeaRobusta|canephora|robusta|rob';

my $run_orthogroups = 1;
my $run_unassigned  = 1;
my $help            = 0;

GetOptions(
  'familywise=s'      => \$familywise,
  'unassigned=s'      => \$unassigned,
  'outdir=s'          => \$outdir,

  'family-cols=s'     => \$family_cols,

  'arabica-col=s'     => \$arabica_col_pat,
  'eug-col=s'         => \$eug_col_pat,
  'rob-col=s'         => \$rob_col_pat,

  'orthogroups!'      => \$run_orthogroups,
  'unassigned-run!'   => \$run_unassigned,

  'help!'             => \$help,
) or die usage();

if ($help) { print usage(); exit 0; }

make_path($outdir) unless -d $outdir;

my @fam_cols = split(/\s*,\s*/, $family_cols);
die "[ERROR] --family-cols must list exactly 6 column names\n" unless @fam_cols == 6;

if ($run_orthogroups) {
  run_unique_orthogroups($familywise, $outdir, \@fam_cols);
}
if ($run_unassigned) {
  run_unique_unassigned($unassigned, $outdir, $arabica_col_pat, $eug_col_pat, $rob_col_pat);
}

print "\n[DONE] Outputs written to: $outdir\n";
exit 0;

# ------------------------------------------------------------
# (A) Unique orthogroups from familywise count table
# ------------------------------------------------------------
sub run_unique_orthogroups {
  my ($in, $outdir, $cols6) = @_;

  open(my $IN, "<", $in) or die "[ERROR] Cannot open familywise table: $in ($!)\n";

  # map column position (1..6) to letters a..f (kept from your original script)
  my %tbl = ( 1=>'a', 2=>'b', 3=>'c', 4=>'d', 5=>'e', 6=>'f' );

  my %ortho_tbl;  # orthoID => [dist_string, total]
  my $header = <$IN>;
  die "[ERROR] Empty familywise file: $in\n" unless defined $header;

  while (my $line = <$IN>) {
    chomp $line;
    next unless length $line;
    my @ele = split(/\t/, $line, -1);
    next if $ele[0] =~ /Orthogroup/i;

    # Expect: Orthogroup + 6 numeric columns
    if (@ele < 7) {
      warn "[WARN] Skipping short line (need 7 cols): $line\n";
      next;
    }

    my $id = $ele[0];
    my $total = sum(map { ($_ =~ /^\d+$/) ? $_ : 0 } @ele[1..6]);

    my $dist = join('', map {
      my $i = $_;
      ($ele[$i] && $ele[$i] =~ /^\d+$/ && $ele[$i] > 0) ? $tbl{$i} : ()
    } (1..6));

    $ortho_tbl{$id} = [$dist, $total];
  }
  close $IN;

  # coffee-only distributions (d/e/f correspond to arabica/eug/robusta in your table order)
  my %unique_groups = (
    'de'  => 'ara.eug',
    'def' => 'ara.eug.rob',
    'd'   => 'ara',
    'df'  => 'ara.rob',
    'f'   => 'rob',
    'ef'  => 'eug.rob',
    'e'   => 'eug',
  );

  foreach my $node (sort keys %unique_groups) {
    my $species = $unique_groups{$node};
    my $out = "$outdir/$species" . "_unique_ortho.tsv";
    open(my $OUT, ">", $out) or die "[ERROR] Cannot write: $out ($!)\n";
    foreach my $ortho (sort keys %ortho_tbl) {
      my $distribution = $ortho_tbl{$ortho}->[0];
      if ($distribution eq $node) {
        print $OUT "$ortho\t$species\n";
      }
    }
    close $OUT;
  }

  print "[OK] Unique orthogroups written (ara/eug/rob combinations)\n";
}

# ------------------------------------------------------------
# (B) Unique unassigned genes from Orthogroups_UnassignedGenes.tsv
# ------------------------------------------------------------
sub run_unique_unassigned {
  my ($in, $outdir, $arab_pat, $eug_pat, $rob_pat) = @_;

  open(my $IN, "<", $in) or die "[ERROR] Cannot open unassigned table: $in ($!)\n";
  my $hdr = <$IN>;
  die "[ERROR] Empty unassigned file: $in\n" unless defined $hdr;
  chomp $hdr;

  my @h = split(/\t/, $hdr, -1);

  # Find coffee columns by header pattern (case-insensitive)
  my $arab_i = find_col_index(\@h, $arab_pat);
  my $eug_i  = find_col_index(\@h, $eug_pat);
  my $rob_i  = find_col_index(\@h, $rob_pat);

  die "[ERROR] Could not find arabica column with pattern: $arab_pat\n" unless defined $arab_i;
  die "[ERROR] Could not find eugenioides column with pattern: $eug_pat\n" unless defined $eug_i;
  die "[ERROR] Could not find robusta/canephora column with pattern: $rob_pat\n" unless defined $rob_i;

  my %genes_arab;
  my %genes_eug;
  my %genes_rob;

  while (my $line = <$IN>) {
    chomp $line;
    next unless length $line;
    my @ele = split(/\t/, $line, -1);
    next if $ele[0] =~ /Orthogroup/i;

    # cells may contain:
    #   - empty
    #   - single gene ID
    #   - multiple gene IDs separated by commas and/or spaces
    add_genes_from_cell(\%genes_arab, $ele[$arab_i] // '');
    add_genes_from_cell(\%genes_eug,  $ele[$eug_i]  // '');
    add_genes_from_cell(\%genes_rob,  $ele[$rob_i]  // '');
  }
  close $IN;

  # Write outputs: one gene per line (tsv: gene_id)
  write_gene_list("$outdir/arabica_unique_unassigned.tsv", \%genes_arab);
  write_gene_list("$outdir/eug_unique_unassigned.tsv",     \%genes_eug);
  write_gene_list("$outdir/robusta_unique_unassigned.tsv", \%genes_rob);

  print "[OK] Unique unassigned genes written (arabica/eug/robusta)\n";
}

sub find_col_index {
  my ($headers, $pattern) = @_;
  my $re = qr/$pattern/i;
  for (my $i=0; $i<@$headers; $i++) {
    return $i if defined $headers->[$i] && $headers->[$i] =~ $re;
  }
  return undef;
}

sub add_genes_from_cell {
  my ($set, $cell) = @_;
  return unless defined $cell;
  $cell =~ s/^\s+|\s+$//g;
  return unless length $cell;

  # split on commas/semicolons/whitespace
  my @ids = split(/[,\s;]+/, $cell);
  for my $id (@ids) {
    next unless defined $id;
    $id =~ s/^\s+|\s+$//g;
    next unless length $id;
    $set->{$id} = 1;
  }
}

sub write_gene_list {
  my ($path, $set) = @_;
  open(my $OUT, ">", $path) or die "[ERROR] Cannot write: $path ($!)\n";
  foreach my $g (sort keys %$set) {
    print $OUT "$g\n";
  }
  close $OUT;
}

sub usage {
  return <<'USAGE';
Usage:
  perl coffee_unique_sets.pl [options]

Inputs:
  --familywise FILE     Familywise count TSV (default: ../1_orthologs_familyWise/0_each_ortho_familywiseCount.tsv)
  --unassigned FILE     Orthogroups_UnassignedGenes.tsv (default: ../Orthogroups_UnassignedGenes.tsv)

Output:
  --outdir DIR          Output directory (default: .)

Familywise config:
  --family-cols STR     Six column names (comma-separated) after Orthogroup
                        default: laminales,solnales,olden,arabica,eug,robusta

Unassigned config (header matching patterns; case-insensitive regex):
  --arabica-col STR     default: CoffeaArabica|arabica
  --eug-col STR         default: CoffeaEugenioides|eugenioides|eug
  --rob-col STR         default: CoffeaRobusta|canephora|robusta|rob

Run toggles:
  --orthogroups / --no-orthogroups        Run/skip unique orthogroup extraction (default: run)
  --unassigned-run / --no-unassigned-run  Run/skip unique unassigned extraction (default: run)

Examples:
  perl coffee_unique_sets.pl --outdir ./out

  perl coffee_unique_sets.pl \
    --familywise ../1_orthologs_familyWise/0_each_ortho_familywiseCount.tsv \
    --unassigned ../Orthogroups_UnassignedGenes.tsv \
    --outdir ./out

  # If headers differ in unassigned file:
  perl coffee_unique_sets.pl --unassigned Orthogroups_UnassignedGenes.tsv \
    --arabica-col 'CoffeaArabica_LP_NoPartial' \
    --eug-col 'CoffeaEugenioides_LP_NoPartial' \
    --rob-col 'CoffeaRobusta_LP_NoPartial'
USAGE
}
