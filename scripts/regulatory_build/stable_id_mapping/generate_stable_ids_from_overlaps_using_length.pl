use strict;
use Data::Dumper;
use Getopt::Long;

my $all_overlaps;
my $old_regulatory_features;
my $new_regulatory_features;
my $outfile;
my $stable_id_prefix;

GetOptions (
'all_overlaps=s'             => \$all_overlaps,
'old_regulatory_features=s'  => \$old_regulatory_features,
'new_regulatory_features=s'  => \$new_regulatory_features,
'outfile=s'                  => \$outfile,
'stable_id_prefix=s'         => \$stable_id_prefix,
);

## Identify the maximum stable id seen in the "old" regulatory feature dataset
open my $old_regulatory_features_fh, '<', $old_regulatory_features;
my $max_seen_stable_id_number = find_max_seen_stable_id( $old_regulatory_features_fh );
$old_regulatory_features_fh->close();

## Retain overlaps between the same regulatory feature types and
## for each one of those only the one with the longest overlap
open my $all_overlaps_fh, '<', $all_overlaps;
my $stable_id_mappings_ref_hash = filter_overlaps( $all_overlaps_fh );
$all_overlaps_fh->close();

open my $new_regulatory_features_fh, "<", $new_regulatory_features;

(
my $stable_id_hash,
my $num_mapped_stable_ids,
my $num_new_stable_ids
)= generate_stable_id_assignments(
$new_regulatory_features_fh,
$stable_id_prefix,
$stable_id_mappings_ref_hash,
$max_seen_stable_id_number
);

$new_regulatory_features_fh->close();
print "Writing stable id assignments to $outfile\n";

open my $out_fh, ">" . $outfile;
foreach my $regulatory_feature_id (keys %$stable_id_hash) {
  
  my $stable_id = $stable_id_hash->{$regulatory_feature_id};
  $out_fh->print(
  join "\t", $regulatory_feature_id, $stable_id
  );
  $out_fh->print("\n");
}
$out_fh->close();

print "$num_mapped_stable_ids stable ids were mapped.\n";
print "$num_new_stable_ids stable ids were newly assigned.\n";


sub find_max_seen_stable_id {
  
  my $bedfile_fh = shift;
  
  my $max_seen_stable_id_number;
  
  while (my $current_bed_file_line = <$bedfile_fh>) {
    # 18	76429380	76430144	Open chromatin	00000105157	535842
    chomp $current_bed_file_line;
    my ($chr, $start, $end, $regulatory_feature_type, $old_regulatory_feature_stable_id, $old_regulatory_feature_db_id) = split "\t", $current_bed_file_line;
    
    # Look for maximum ID number among the old features
    if ($max_seen_stable_id_number < $old_regulatory_feature_stable_id) {
      $max_seen_stable_id_number = $old_regulatory_feature_stable_id;
    }
  }
  
  print "The maximum stable id seen is: $max_seen_stable_id_number\n";
  return( $max_seen_stable_id_number );
}

sub filter_overlaps {
  
  my $all_overlaps_fh = shift;
  
  my $same_type_overlaps = retain_same_feature_type_overlaps($all_overlaps_fh);
 
  my %old_pref = ();
  my %new_pref = ();

  my $no_overlaps_retained=0;
  my $no_overlaps_dismissed=0;
  
  foreach my $overlap ( sort { $b->[12] <=> $a->[12] } @{$same_type_overlaps} ) {
    my $old_regulatory_feature_stable_id  = $overlap->[4];
    my $new_regulatory_feature_type_db_id = $overlap->[11];
    
    if( ! (exists $old_pref{$old_regulatory_feature_stable_id} || exists $new_pref{$new_regulatory_feature_type_db_id} ) ) {
      $old_pref{$old_regulatory_feature_stable_id} = $new_regulatory_feature_type_db_id;
      $new_pref{$new_regulatory_feature_type_db_id} = $old_regulatory_feature_stable_id;
      
      $no_overlaps_retained += 1;
      
    } else {
      
      $no_overlaps_dismissed += 1;
      #print "Either the $old_regulatory_feature_stable_id or $new_regulatory_feature_type_db_id have already been used\n";
    }
  }
  
  print "$no_overlaps_retained overlap were retained.\n";
  print "$no_overlaps_dismissed overlaps were discarded.\n";
  
  return ( \%new_pref );
}

sub retain_same_feature_type_overlaps {
  my $bedfile_fh = shift;
  
  my @same_type_overlaps;

  while (my $current_bed_file_line = <$bedfile_fh>) {
    #1	13201	13800	Enhancer	00000000001	623456	1	13371	13724	Open chromatin	00000341931	160079	353
    
    chomp $current_bed_file_line;
    my @bed_file_fields = split "\t", $current_bed_file_line;
    
    my $old_regulatory_feature_type = $bed_file_fields[3];
    my $new_regulatory_feature_type = $bed_file_fields[9];
    
    # Feature types must be identical for transferring the stable id.
    if ($old_regulatory_feature_type eq $new_regulatory_feature_type) {
      #print "Compatible $current_bed_file_line\n";
      push @same_type_overlaps, \@bed_file_fields;
    } #else {
      #print "Not compatible $current_bed_file_line\n";
    #}
  }
  
  return( \@same_type_overlaps );
}

sub generate_stable_id_assignments {
  
  my $new_regulatory_features_fh  = shift;
  my $stable_id_prefix            = shift;
  my $stable_id_mappings_ref_hash = shift;
  my $max_seen_stable_id_number   = shift;
  
  my %stable_id_mappings_hash = %$stable_id_mappings_ref_hash;
  
  my @no_stable_id_mappings = keys %$stable_id_mappings_ref_hash;
  print "Recovered ". @no_stable_id_mappings." mappings!!\n";
  
  my %stable_id_hash = ();
  my $next_free_id = $max_seen_stable_id_number + 1;
  
  my $num_mapped_stable_ids = 0;
  my $num_new_stable_ids    = 0;
  
  while (my $line = <$new_regulatory_features_fh>) {
    #15	102118695	102119230	TF binding site	00000368862	1
    chomp $line;
    my @bed_file_fields = split "\t", $line;
    
    my $updated_stable_id = undef;
    my $new_regulatory_feature_type_db_id = $bed_file_fields[5];
    
    if ( exists $stable_id_mappings_hash{$new_regulatory_feature_type_db_id} ) {
      $updated_stable_id = $stable_id_mappings_hash{$new_regulatory_feature_type_db_id};
      $num_mapped_stable_ids++;
    } else {
      $updated_stable_id = $next_free_id;
      $next_free_id += 1;
      $num_new_stable_ids++;
    }
    
    # Creating stable id string, composed of prefix + 11 digit integer, front padded with 0s
    $stable_id_hash{$new_regulatory_feature_type_db_id} = $stable_id_prefix . sprintf("%011d", $updated_stable_id);
    
  }
  return \%stable_id_hash, $num_mapped_stable_ids, $num_new_stable_ids;
}





