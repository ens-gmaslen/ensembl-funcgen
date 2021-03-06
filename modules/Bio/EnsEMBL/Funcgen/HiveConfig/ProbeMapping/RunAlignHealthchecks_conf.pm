package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::RunAlignHealthchecks_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    my $max_probe_features_to_check = 10000;
    
    return [
      {
          -logic_name  => 'start_align_healthchecks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'check_duplicate_probe_features'
          },
      },
      {   -logic_name        => 'check_duplicate_probe_features',
          -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => '
                delete_redundant_probe_features.pl \
                  --registry #reg_conf# \
                  --species #species# \
                  --only_test 1
              ',
          },
          -flow_into => {
              MAIN => 'check_gapped_probe_features_from_transcript_matches'
          },
      },
      {
          -logic_name  => 'check_gapped_probe_features_from_transcript_matches',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 'check_probe_feature_sequences.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --logic_name ProbeAlign_transcript '
                . ' --check_probe_features_with_nontrivial_cigar_lines 1'
                . ' --max_check ' . $max_probe_features_to_check
          },
          -flow_into => {
              MAIN => 'check_ungapped_probe_features_from_transcript_matches',
          },
        },
      {
          -logic_name  => 'check_ungapped_probe_features_from_transcript_matches',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 'check_probe_feature_sequences.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --logic_name ProbeAlign_transcript '
                . ' --max_check ' . $max_probe_features_to_check
          },
          -flow_into => {
              MAIN => 'check_probe_feature_sequences_from_genomic_matches',
          },
        },
      {
          -logic_name  => 'check_probe_feature_sequences_from_genomic_matches',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 'check_probe_feature_sequences.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --logic_name ProbeAlign_genomic '
                . ' --max_check ' . $max_probe_features_to_check
          },
        },
    ];
}

1;
