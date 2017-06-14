use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
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

open my $old_regulatory_features_fh, '<', $old_regulatory_features;

my $max_seen_stable_id_number = find_max_seen_stable_id($old_regulatory_features_fh);
$old_regulatory_features_fh->close();

open my $all_overlaps_fh, '<', $all_overlaps;

my $overlaps = filter_for_compatible_overlaps($all_overlaps_fh);
$all_overlaps_fh->close();

(
  my $stable_id_hash,
  my $num_mapped_stable_ids,
  my $num_new_stable_ids  
)= generate_stable_id_assignments(
  $stable_id_prefix, 
  $overlaps, 
  $max_seen_stable_id_number
);

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

sub generate_stable_id_assignments {

  my $stable_id_prefix          = shift;
  my $overlaps                  = shift;
  my $max_seen_stable_id_number = shift;

  # Go through overlaps in order of increasing overlap length. This means that you should always
  # overwrite an overlap with a later one.
  # 
  # The overlap length is in column 12, so the index is 11.
  sub cmp_overlaps {
    return $a->[11] <=> $b->[11];
  }
  
  my %old_pref = ();
  my %new_pref = ();
  foreach my $entry (sort cmp_overlaps @{$overlaps}) {
    $old_pref{$entry->[4]} = $entry->[10];
    $new_pref{$entry->[10]} = $entry->[4];
  }

  my %stable_id_hash = ();
  my $next_free_id = $max_seen_stable_id_number + 1;
  open my $in, "<", $new_regulatory_features;

  my $num_mapped_stable_ids = 0;
  my $num_new_stable_ids    = 0;

  while (my $line = <$in>) {
    chomp $line;
    my @items = split "\t", $line;
   
    #print "Compatible $line\n";
 
    # This is the name of the feature in the bed file. e.g.: tfbs_2, ctcf_1, 
    # dnase_2, dnase_3, tfbs_3, tfbs_4, proximal_1, proximal_2, tss_1, 
    # tfbs_6, tfbs_7, ctcf_2, ctcf_3, tfbs_8 etc.
    #
    my $new_stable_id_number = $items[4];

    my $stable_id = undef;

    if (exists $new_pref{$new_stable_id_number} && $old_pref{$new_pref{$new_stable_id_number}} eq $new_stable_id_number) {
      $stable_id = $new_pref{$new_stable_id_number};
      $num_mapped_stable_ids++;
    } else {
      $stable_id = $next_free_id;
      $next_free_id += 1;
      $num_new_stable_ids++;
    }

    # Creating stable id string, composed of prefix + 11 digit integer, front padded with 0s
    # 
    my $regulatory_feature_id = $items[5];
    $stable_id_hash{$regulatory_feature_id} = $stable_id_prefix . sprintf("%011d", $stable_id);

  }
  close $in;
  return \%stable_id_hash, $num_mapped_stable_ids, $num_new_stable_ids;
}

sub filter_for_compatible_overlaps {

  my $all_overlaps_fh = shift;
 
  my @overlaps;
  my %overlaps_length;

  my $all_overlaps = read_all_overlaps( $all_overlaps_fh );

  foreach my $entry (sort { $b->[12] <=> $a->[12] } @{$all_overlaps}) {
    my $overlap_length = $entry->[12];
    my $regulatory_feature_type_old  = $entry->[4];

    if( ! exists ($overlaps_length{$regulatory_feature_type_old}) ) {
      push @overlaps, $entry;
      $overlaps_length{$regulatory_feature_type_old} = $overlap_length;
      #print "Looking at entry: @{$entry}\n";
    } elsif ( $overlaps_length{$regulatory_feature_type_old} < $overlap_length) {
      print "ERROR: The new entry has a longer overlap than the one already reported.\n@{$entry}\n";
    }
  }
  return ( \@overlaps );
}

sub find_max_seen_stable_id {
    
  my $bedfile_fh = shift;
    
  my $max_seen_stable_id_number;
    
  while (my $current_bed_file_line = <$bedfile_fh>) {
      
    chomp $current_bed_file_line;
    my @bed_file_fields = split "\t", $current_bed_file_line;
      
    my $regulatory_feature_old_stable_id_number = $bed_file_fields[4];
    
      # Look for maximum ID number among the old features
    if ($max_seen_stable_id_number < $regulatory_feature_old_stable_id_number) {
        $max_seen_stable_id_number = $regulatory_feature_old_stable_id_number;
    }
  }

  print "The maximum stable id seen is: $max_seen_stable_id_number\n";
  return( $max_seen_stable_id_number );
}

sub read_all_overlaps {
  my $bedfile_fh = shift;
    
  my @all_overlaps;

  while (my $current_bed_file_line = <$bedfile_fh>) {
        
    chomp $current_bed_file_line;
    my @bed_file_fields = split "\t", $current_bed_file_line;
    
    my $regulatory_feature_type_old  = $bed_file_fields[3];
    my $regulatory_feature_type_new  = $bed_file_fields[9];
      
    if ($regulatory_feature_type_old eq $regulatory_feature_type_new) {
      # Feature types must be identical for transferring the stable id.
      
      #print "Compatible $current_bed_file_line\n";
      push @all_overlaps, \@bed_file_fields;
    } #else {
      #print "Not compatible $current_bed_file_line\n";
    #}
  }
  return( \@all_overlaps );
}
