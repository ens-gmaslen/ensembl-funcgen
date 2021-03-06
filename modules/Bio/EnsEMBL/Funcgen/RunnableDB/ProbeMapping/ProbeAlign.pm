# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::ProbeMapping::ProbeAlign;




=head1 SYNOPSIS

  my $affy = Bio::EnsEMBL::Analysis::RunnableDB::ProbeMapping::ProbeAlign->new
    (-db         => $refdb,
     -analysis   => $analysis_obj,
     -database   => $EST_GENOMIC,
     -query_seqs => \@sequences,
     );

  $affy->fetch_input();
  $affy->run();
  $affy->write_output(); #writes to DB

=head1  DESCRIPTION

This object maps probes to a genomic and transcript sequence, writing the results as 
Bio::EnsEMBL::Funcgen::ProbeFeatures. UnmappedObjects are also written for those probes
which either do not map at all or exceed the maximum mapping threshold defined by HIT_SATURATION_LEVEL.
You must FIRST have created and imported all the necessary Arrays and Probe objects using the ImportArrays module.

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>

=head1 METHODS

=cut

#To do
#1 Enable rollback based on the probe IDs for this chunk
#2 Delay Unmapped object write? 
#This won't be necessary if we enable rollback and will reduce mem usage
#If we can write features immediately too.
#3 Report Unmapped object count in output
#4 Re/move all unnecessary warns

package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
# use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ExonerateProbe;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBEntry;

#use Bio::EnsEMBL::Analysis::Config::ProbeAlign;

#use vars qw(@ISA);
#@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB);
use base ('Bio::EnsEMBL::Hive::Process');

############################################################


#TO DO
# 1. Account for 5'/3' hanging alignment mismatches
#    These may be a sequence match or mismatch, but are currently being reported as an alignment 'm'
#    Which implies a sequence mis-match. This is simple as we don't generally have long overhangs
#    due to restrictive mismatch rules in ExonerateProbe
#    But would get more complicated (i.e. performing alignment) with longer overhangs
#    Can we not just report the shorter alignment?  This will cause problems when trying 
#    to identify which bp of the probe a particular mismatch is on i.e.13th for AFFY MM probe?
#    Can we substitute M for E(nd) or U(nknown)? This would have to be done in ExonerateProbe
#    Is potentially a source for problems if a probe were to overhand the 3' end and this were not
#    caught by the genomic mapping due to an intron being present at the very end of the transcript.
#    e.g. ENSMUST00000097969
#    To capture this we must perform 5'/3' extension on the alignment against the genome!
#    Would need to get the query sequence from ExonerateProbe or directly from the fasta
#    Change to U(nknown not exonerate standard nomenclature) or N(on-equivalenced region), the
#    later is normally used for large regions on low conservation in protein alignments e.g. loops
#    THIS IS NOW SET TO U IN set_probe_and_slice
#
# 2. Delay writes of unmapped objects so we don't write anything unless -w is specified when testing
#    Also prevents having to rollback if job failed hlafway through filter
#
# 3. Enable rollback by chunk.  This should skip the run

sub new {
  my ( $class, @args ) = @_;
  my $self = bless {},$class;  
  return $self;
}

sub probe_align_helper {
  my ($self, $probe_align_helper) = @_;
  if ($probe_align_helper) {
    confess('Type error') unless ($probe_align_helper->isa('Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlignHelper'));
    $self->{'probe_align_helper'} = $probe_align_helper;
  }
  return $self->{'probe_align_helper'};
}

# sub mapping_target_type {
#   my ($self, $mapping_target_type) = @_;
#   if ($mapping_target_type) {
#     $self->{'mapping_target_type'} = $mapping_target_type;
#   }
#   return $self->{'mapping_target_type'};
# }

sub fetch_input {

  my $self = shift;
  
  $self->IIDREGEXP('(\d+):(\d+)');
  
  
  # We need 101 hits, so we can test later, if a probe sequencehad more than 100 hits.
  $self->OPTIONS(' --bestn 101 --dnahspthreshold 116 --fsmmemory 256 --dnawordlen 12 --dnawordlimit 0 ');
  $self->HIT_SATURATION_LEVEL(100);
  $self->MAX_MISMATCHES(0);
  $self->QUERYTYPE('dna');
  
  my $species                 = $self->param('species');
  my $funcgen_adaptor         = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

  my $out_db = 
  
  $self->OUTDB($funcgen_adaptor);
  $self->QUERYSEQS($self->param('QUERYSEQS'));
  
  use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils qw( create_dna_db_params_from_funcgen_hash );    
  my %dna_db_hash = create_dna_db_params_from_funcgen_hash($self->param('OUTDB'));
  
  $self->DNADB(\%dna_db_hash);
  #$self->mapping_target_type($self->param('mapping_target_type'));
  #$self->mapping_target_type('transcript');
  
  my $tracking_dba_hash = $self->param('tracking_dba_hash');
  use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils qw( create_funcgen_adaptor );
  my $funcgen_adaptor = create_funcgen_adaptor($tracking_dba_hash);
  
  $funcgen_adaptor->dbc->reconnect_when_lost(1);
  $funcgen_adaptor->dnadb->dbc->reconnect_when_lost(1);  
  $funcgen_adaptor->dbc->disconnect_when_inactive(1);
  $funcgen_adaptor->dnadb->dbc->disconnect_when_inactive(1);
  
  $self->dbc->disconnect_when_inactive(1);
  
  my $analysis_adaptor = $funcgen_adaptor->get_AnalysisAdaptor;
  
  use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlignHelper;    
  my $probe_align_helper = Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlignHelper->new(
    -funcgen_adaptor => $funcgen_adaptor,
  );
  
  $self->probe_align_helper($probe_align_helper);
  
  #Set up and check DBs before we do the mapping
#   $self->outdb->dbc->db_handle;
#   $self->outdb->dnadb->dbc->db_handle;

  #Set session server mode to prevent truncated cigar_lines
  $self->outdb->dbc->disconnect_when_inactive(1);
  $self->outdb->dbc->db_handle->do('set SESSION sql_mode="STRICT_ALL_TABLES"');

  $self->arrays({});
  $self->probes({});
  
  my $mapping_type = $self->param('mapping_type');
  
  my $schema_build = $self->outdb->_get_schema_build($self->outdb->dnadb);
  my $species = $self->outdb->species;
  #throw('Must provide a -species to the OUTDB in Bio::EnsEMBL::Analysis::Config::ProbeAlign') if ! $species;
	 
  my ($db_name, $display_name);

  #Check imported status of arrays
  #($array_class = $logic) =~ s/_Probe.*Align$//;

  #The IMPORTED status is at the ArrayChip level not the Array level!
  #$self->outdb->get_ArrayAdaptor->check_status_by_class('IMPORTED', $array_class);


  if($mapping_type eq 'transcript'){
	 $self->{'mapping_type'} = 'transcript';
	 
	 #We need to check for ensembl_Core_Transcript external db here
	 #This should also be linked to the release
	 #Should we store here? May over-write each other?
	 #No problem if it fails once as the pipeline will resubmit the job
	 
	 #Will this cause problems with maintaining a master external_dbs file
	 #As this may will grow *3 for every release.
	 
	 #But once we allow GENE, TRANSCRIPT, PROTEIN in xref.info_type
	 #Then this would only be one extra db per release
	 
	 #Would be DBEntry::ADaptor->get_all_external_dbs_by_name

	 #This is tricky as we may be aligning against an old DB
	 #i.e. if there was a new array but not a new build
	 #Current DBEntry queries are naive to external_db version

	 $db_name      = $species.'_core_Transcript';
	 $display_name = 'EnsemblTranscript';
  }
  else{
	 $self->{'mapping_type'} = 'genomic';
	 $db_name      = $species.'_core_Genome';
	 $display_name = 'EnsemblGenome';	 

	 #This is really just a place holder for UnmappedObjects against the entire genome
	 #Should we redo this as EnsemblSpecies and actually have xref entries for each species name?
  }


  my $sql = "select external_db_id from external_db where db_name='$db_name' and db_release='$schema_build'";
  my ($extdb_id) = $self->outdb->dbc->db_handle->selectrow_array($sql);
	
 
  #This is causing redundant edb entries as parallel job store the same entry, which is not constrained by unique key.
  #Need to store this before hand and fail here

  if(! $extdb_id){
	#status is dubious here as this should really be on object_xref
	my $insert_sql = 'INSERT into external_db(db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type)'.
	  " values('$db_name', '$schema_build', 'KNOWNXREF', 1, 5, '$display_name', 'MISC')";
	throw("Failed to fetch external_db_id for $db_name $schema_build\nPlease add this record using the following sql:\n$insert_sql");


  }

  $self->{'_external_db_id'} = $extdb_id;

  ##########################################
  # set up the target (genome)
  ##########################################

  $self->TARGETSEQS($self->param('TARGETSEQS'));
  
  my $target = $self->TARGETSEQS;

  if ( -e $target ){
    if(-d $target ) {
      print ("Target $target is a directory of files\n");
    }elsif (-s $target){
      print ("Target $target is a whole-genome file\n");
    }else{
      throw("'$target' isn't a non-empty file or a directory");
    }
  } else {
    throw("'$target' target seqs could not be found");
  }

  ##########################################
  # set up the query probe seq
  ##########################################

  my ($query_file, $chunk_number, $chunk_total);

  my $query = $self->QUERYSEQS;
  
  if ( -e $query and -s $query ) {

    # query seqs is a single file; input id will correspond to a chunk number
    $query_file = $query;
    my $iid_regexp = $self->IIDREGEXP;

    if (not defined $iid_regexp){
      throw("You must define IIDREGEXP in config to enable inference of chunk number and total from your single fasta file" )
    }

    #( $chunk_number, $chunk_total ) = $self->input_id =~ /$iid_regexp/;
#     $chunk_number = $self->param('chunk_number');
#     $chunk_total  = $self->param('chunk_total');
#     
#     if(!$chunk_number || !$chunk_total){
#       throw "I can't make sense of your input id  using the IIDREGEXP in the config!\n";
#     }

    #store this for reference later
    $self->query_file($query_file);

  } else {

    throw("'$query'  must refer to a single fasta file with all probe sequences referenced by affy_probe_id\n");

  }

  ##########################################
  # setup the runnable
  ##########################################

  #my %parameters = %{ $self->parameters_hash };
  my %parameters;
 
  #Parameter hash is picked up from the DB analysis table
  #else the config -options value is used
  

  if (not exists $parameters{-options}) {
    if ($self->OPTIONS) {
      $parameters{-options} = $self->OPTIONS;
    } else {
	  throw('You have not options stored in the DB or accessible from the analysis config hash for '.$self->analysis->logic_name);
      #$parameters{-options} = "";
    }
  }
  
  #print STDERR "PROGRAM FILE: ".$self->analysis->program_file."\n";

  my $runnable = Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ExonerateProbe->new
	(
# 	 -program            => $self->analysis->program_file,
# 	 -analysis           => $self->analysis,
	 -program            => 'exonerate',
	 -analysis           => 'Where is this used',
	 -target_file        => $target,
	 -query_type         => $self->QUERYTYPE,#dna
	 -query_file         => $query_file,
	 -query_chunk_number => $chunk_number,
	 -query_chunk_total  => $chunk_total,
	 -max_mismatches     => $self->MAX_MISMATCHES,
	 -mapping_type       => $self->mapping_type,
	 #-filter_method      => $self->FILTER_METHOD,

	 probe_align_helper =>  $self->probe_align_helper,
	 
	 %parameters,
  );
  
  $self->runnable($runnable);

}

sub runnable{
  my ($self, $runnable) = @_;
  if(!$self->{'runnable'}){
    $self->{'runnable'} = [];
  }
  if($runnable){
    push(@{$self->{'runnable'}}, $runnable);
  }
  return $self->{'runnable'};
}

sub mapping_type{
  return $_[0]->{'mapping_type'};
}

############################################################
#Overrides RunnableDB::run

sub run {
  my ($self) = @_;

  throw("Can't run - no runnable objects") unless ( $self->runnable );

  my $runnable = @{ $self->runnable }[0];
  $runnable->run;


  #Can we rollback individual chunk output here
  #so we never duplicate output
  #If job fails we currently don't know whether it failed before storing any unmapped objects
  #Can we do this optionally by setting an env var?
  #We would set this in SubmitAlign via a parameter
  #We need to set this in the config to maintain decoupling between the environment and the RunnableDB
  #Could we record a status in the eFG DB to track this?
  #Don't need both?
  #does not fit with efg status as is not a table, but could be hijacked quite easily?
  #Let's just set this manually. You will have to reset the jobs, so we could add a warning here to use the appropriate flag?
  


  
  #These must be used in the rule_manager.pl
  #Which is where write_output must also be called ffrom
  #$self->output($features);#Moved this after filter/set_probe so we get the write output reported
  $self->filter_features($runnable->output);
}

############################################################
#overrides RunnableDB::write_output

sub write_output {
  my ( $self, @output ) = @_;

  my $outdb           = $self->outdb;
  my $probe_adaptor   = $outdb->get_ProbeAdaptor;
  my $feature_adaptor = $outdb->get_ProbeFeatureAdaptor;
  my $dbe_adaptor     = $outdb->get_DBEntryAdaptor;
  
  #Add analysis, slices to features, and make
  #sure they're pointing at the persistent array / probe instances
  #instead of the fake arrays & probes they were created with
  #Can we not do this in filter as we are not writing here?
  #my $features = $self->set_probe_and_slice($self->features);

  #Now set correct output
  #This is only used to report the number of features
  my $features = $self->features;
  $self->output($features);
  
  print 'Writing '.scalar(@{$features})." ProbeFeatures\n";
  
  $outdb->dbc->disconnect_when_inactive(0);

  foreach my $feature_xref(@{$features}){
  
    my ($feature, $xref) = @$feature_xref;

    if(! defined $feature->slice()){
      warn "Skipping feature storage for Probe, no slice available for ".$feature->seqname();
      next;
    }

#     my @probe_analysis_pair = $self->probe_align_helper->fetch_probe_and_analysis_for_probe_feature_linking(
#       $feature->probe_id,
#       'genomic'
#     );
#     foreach my $current_probe_analysis_pair (@probe_analysis_pair) {
#     
#       my $probe    = $current_probe_analysis_pair->{probe};
#       my $analysis = $current_probe_analysis_pair->{analysis};
#       
#       foreach my $current_analysis (@$analysis) {
#       
# 	use Bio::EnsEMBL::Funcgen::ProbeFeature;
# 	
# 	# Create a new feature from scratch to make sure the store method
# 	# can't set database ids which would mess this up.
# 	#
# 	my $feature_per_probe = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
# 		-PROBE         => $probe,
# 		-MISMATCHCOUNT => $feature->mismatchcount,
# 		-SLICE         => $feature->slice,
# 		-START         => $feature->start,
# 		-END           => $feature->end,
# 		-STRAND        => $feature->strand,
# 		-ANALYSIS      => $current_analysis,
# 		-CIGAR_STRING  => $feature->cigar_string,
# 	);

	eval { 
	  $feature_adaptor->store($feature)
	};
        if ($@) {

          # Don't terminate for duplicate features. This can happen, if a job 
          # is restarted by hive.
          #
          if ($@!~/st execute failed: Duplicate entry/) {
            $self->throw('Unable to store ProbeFeature for probe '.$feature->probe_id." on slice:\t".$feature->slice->name."\n$@");
          }
        }
	if($xref){
	  #Remove ignore release flag so we store on the correct schema build
	  eval{ $dbe_adaptor->store($xref, $feature->dbID, 'ProbeFeature') };

	  if ($@) {
		$self->throw('Unable to store ProbeFeature DBEntry for probe '.$feature->dbID." !\n$@");
	  }
	} # of if xref      
#       } # foreach analysis    
#     } # of foreach probe analysis pair
  }
}

############################################################

sub filter_features {
  my ($self, $output) = @_;
  
  my $features = $output->{features};
  my $promiscuous_probes = $output->{promiscuous_probes};
  
  use Data::Dumper;  
  my (%hits_by_probe, @kept_hits);
  
  my $mapping_type = $self->mapping_type;
  my $max_hits     = $self->HIT_SATURATION_LEVEL;
  my $uo_adaptor   = $self->outdb->get_UnmappedObjectAdaptor;
  
  # In ExonerateProbe, any probe that makes more than 100 hits is classifies as a promiscuous probe.
  #
  foreach my $promiscuous_probe (@$promiscuous_probes) {
  
    my $probe    = $promiscuous_probe->{probe};
    my $analysis = $promiscuous_probe->{analysis};
    
    foreach my $current_analysis (@$analysis) {
    
      $uo_adaptor->store(
	Bio::EnsEMBL::UnmappedObject->new (
	    -type                => 'array_mapping',
	    -analysis            => $current_analysis,
	    -ensembl_id          => $probe->dbID,
	    -ensembl_object_type => 'Probe',
	    -external_db_id      => $self->{'_external_db_id'},
	    -identifier          => $mapping_type,
	    -summary             => 'Promiscuous probe',
	    -full_desc           => "Probe exceeded maximum allowed number of $mapping_type mappings. Yihaa!"
	)
      );
    }
  }
  
  my $transcript_mapping;
  $transcript_mapping = 1 if $mapping_type eq 'transcript';
  
  foreach my $hit (@$features) {
	push @{$hits_by_probe{$hit->probe_id}}, $hit;
  }

  PROBE_HIT: foreach my $probe_id (keys %hits_by_probe) {
    my @hits = @{$hits_by_probe{$probe_id}};
	my $num_hits = scalar(@hits);
	my $all_hits = $num_hits;	
	
	if ($transcript_mapping) {
	
	  my $maps_promiscuously_to_genome;
	  my @uos = @{$uo_adaptor->fetch_all_by_object_type_id('Probe', $probe_id)};
	  	  
	  UNMAPPED_OBJECT_ENTRY: foreach my $uo(@uos){
		if($uo->summary eq 'Promiscuous probe') {
		  $maps_promiscuously_to_genome = 1;
		  last UNMAPPED_OBJECT_ENTRY;
		}
	  }
	  if ($maps_promiscuously_to_genome) {
	    next PROBE_HIT;
	  }
	}

	if ($num_hits > $max_hits) {
	      $uo_adaptor->store(
		Bio::EnsEMBL::UnmappedObject->new (
		    -type                => 'array_mapping',
		    -analysis            => $hits[0]->analysis,
		    -ensembl_id          => $probe_id,
		    -ensembl_object_type => 'Probe',
		    -external_db_id      => $self->{'_external_db_id'},
		    -identifier          => $mapping_type,
		    -summary             => 'Promiscuous probe',
		    -full_desc           => "Probe exceeded maximum allowed number of $mapping_type mappings(${num_hits}/${max_hits})"
		)
	      );
	      next PROBE_HIT;
	}
	
	if($num_hits == 0) {
	  $uo_adaptor->store(
	    Bio::EnsEMBL::UnmappedObject->new (
	      -type                => 'array_mapping',
	      -analysis            => $hits[0]->analysis,
	      -ensembl_id          => $probe_id,
	      -ensembl_object_type => 'Probe',
	      -external_db_id      => $self->{'_external_db_id'},
	      -identifier          => $mapping_type,
	      -summary             => 'Unmapped probe',
	      -full_desc           => "Probe has no $mapping_type mappings",
	    )
	  );
	  next PROBE_HIT;
	}

	if (
	  ($num_hits > 0) 
	  && ($num_hits <= $max_hits)
	) {
	  push @kept_hits, @hits;
	  next PROBE_HIT;
	}
	# We should never get here.
  }
  $self->set_probe_and_slice(\@kept_hits);
  return \@kept_hits;
}

############################################################

=head2 get_display_name_by_stable_id

  Args [0]   : stable ID from core DB.
  Args [1]   : stable feature type e.g. gene, transcript, translation
  Example    : $self->validate_and_store_feature_types;
  Description: Builds a cache of stable ID to display names.
  Returntype : string - display name
  Exceptions : Throws is type is not valid.
  Caller     : General
  Status     : At risk

=cut

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> display_name
# Return hashref keyed on {$type}{$stable_id}
#Need to update cache if we're doing more than one 'type' at a time
# as it will never get loaded for the new type!


#Replace this with Helper/EFGUtils method

sub get_display_name_by_stable_id{
  my ($self, $stable_id, $type) = @_;

  $type = lc($type);

  if($type !~ /(gene|transcript|translation)/){
	throw("Cannot get display_name for stable_id $stable_id with type $type");
  }
  
  if(! exists $self->{'display_name_cache'}->{$stable_id}){
	($self->{'display_name_cache'}->{$stable_id}) = $self->outdb->dnadb->dbc->db_handle->selectrow_array("SELECT x.display_label FROM $type t, xref x where t.display_xref_id=x.xref_id and t.stable_id='${stable_id}'");
  }

  return $self->{'display_name_cache'}->{$stable_id};
}


sub set_probe_and_slice {
  my ( $self, $features ) = @_;

  print "Setting probe and slice objects\n";
  
  
  $self->{'efg_db'} = undef;

  my $db = $self->outdb;
  
#  $db->RECONNECT_WHEN_CONNECTION_LOST
  
  #$db->dbc->disconnect_when_inactive(1);
  
  my $slice_adaptor = $db->get_SliceAdaptor;
  my $probe_adaptor = $db->get_ProbeAdaptor;
  my $trans_adaptor = $db->dnadb->get_TranscriptAdaptor;
  my $mapping_type  = $self->mapping_type;
  my $gene_adaptor  = $db->dnadb->get_GeneAdaptor;
  #my $dbe_adaptor   = $db->get_DBEntryAdaptor;
  my $analysis      = $self->analysis;
  my $schema_build  = $self->outdb->_get_schema_build($self->outdb->dnadb);
  my $edb_name      = $self->outdb->species.'_core_' . ucfirst($mapping_type);

#   $slice_adaptor->dbc->reconnect(1);
#   $probe_adaptor->dbc->reconnect(1);
#   $trans_adaptor->dbc->reconnect(1);
#   $gene_adaptor ->dbc->reconnect(1);

  my (%slices, $slice_id, $trans_mapper, $align_type, $align_length, @tmp, $gap_length);
  my (@trans_cigar_line, @genomic_blocks, $cigar_line, $gap_start, $block_end);
  my ($block_start, $slice, @gaps, @features, %gene_hits, $gene, @stranded_cigar_line, $gene_sid);
  my ($query_perc, $display_name, $gene_hit_key, %transcript_cache);

  foreach my $feature (@$features) {
    # get the slice based on the seqname stamped on in the runnable
    my $seq_id = $feature->seqname;
	my $probe_id = $feature->probe_id;
	my $load_feature = 1;
	my $xref;
	my ($genomic_start, $genomic_end);

	#Get the slice

	#warn "cigar is ".$feature->cigar_string;

	if($mapping_type eq 'transcript'){

	  if(! exists $transcript_cache{$seq_id}){
		$transcript_cache{$seq_id} = $trans_adaptor->fetch_by_stable_id($seq_id);
	  }

	  print "\n$seq_id\n";
	  $slice_id  =  $transcript_cache{$seq_id}->seq_region_name;

	  if (! exists $slices{$slice_id} ) {
		$slices{$slice_id} = 	$transcript_cache{$seq_id}->slice;
	  }
	}
	else{
	  $slice_id= $seq_id;
	  
	  if ( not exists $slices{$slice_id} ) {
		# assumes genome seqs were named in the Ensembl API Slice naming
		# convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
        $slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
        
        #Temporary hack to get around incorrect header format issue from new sequence_dump
        #$slices{$slice_id} = $slice_adaptor->fetch_by_region(undef, $slice_id);
	  }
	}

    $slice = $slices{$slice_id};
	

	#Now project onto genomic coords if we are on cdna
	if($mapping_type eq 'transcript'){
	  #This will basically need to change the start and stops accordingly
	  #and insert D's the size of introns
	  #reject id no D's found as we will already have this as genomic mapping.

	  #The cigar line should always be displayed with reference to the +ve strand of the
	  #target feature i.e. the -ve strand if the transcript is on the -ve strand
	  #These are local stranded start and ends with respect to transcript strand
	  my $transcript_start  = $feature->start;
	  my $transcript_end    = $feature->end;
	  my $transcript_strand = $transcript_cache{$seq_id}->feature_Slice->strand;
	
	  #warn "\n\n\ntrans($seq_id) start end strand $transcript_start $transcript_end $transcript_strand";


	  $trans_mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript_cache{$seq_id});
	  @genomic_blocks = $trans_mapper->cdna2genomic($transcript_start, $transcript_end);
	

	  #Next feature if this is an ungapped alignment (1 block)
	  #or just representing a flank seq overhang
	  #which will have been caught by the genomic mapping (2 blocks)

	  next if(scalar(grep{$_->isa("Bio::EnsEMBL::Mapper::Coordinate")} @genomic_blocks) < 2);
	  
	  #Coordinate object returned are genomic coords
	  #Introns(ProbeFeature deletions) are only represented by absent Coordindate blocks
	  #5'/3' overhangs(Genomic deletions) are represented by Gap objects
 
	  #We are not accounting for 5'/3' hanging alignment mismatches
	  #Which may actually be a sequence match or mismatch to the genomic sequence!

	  #else alter the start stop values and rebuild the cigarline
	  my $cigar_line = '';
	  @trans_cigar_line = split/:/, $feature->cigar_string;
		 
	  if($transcript_strand == -1){
		#Always report start end and cigar line wrt the +ve strand of the target seq.
		#Even if the match is against the -ve strand. This is what we do in compara.
		#This is always the genome even for transcript mapping, as we a re projecting back.
		#Need to flip genomic blocks around so that we're still dealing with +ve order
		#This is counterintuitive to how the rest of the API deals with stranded features
		#Normally returns lowest start first wrt to +ve strand even if you use a -ve strand slice, no?
		@genomic_blocks      = reverse(@genomic_blocks);
		@stranded_cigar_line = reverse(@trans_cigar_line);

		#This only makes a difference if there are any m's (seq mismatches)
		#print "Transcript is -ve strand so reverse cigar is @stranded_cigar_line\n";


		#This would be an anti-sense transcript hit???

	  }
	  else{
		@stranded_cigar_line = @trans_cigar_line;
	  }


	  #Account for 5'/3' overhanging gaps here?
   
	  my %gap_lengths = (
						 5 => undef,
						 3 => undef,
						);
    #warn "\n\n\ntrans($seq_id) start end strand $transcript_start $transcript_end $transcript_strand";
    #warn "probe if $probe_id with cigar line @stranded_cigar_line";
      

	  foreach my $block(@genomic_blocks){
      #warn $block.' '.$block->start.' '.$block->end;#.' '.$block->strand;
		
      if(! $genomic_start){
        
        if($block->isa('Bio::EnsEMBL::Mapper::Coordinate')){#Set genomic_start
          
          if($gap_lengths{5}){#We have seen a gap
            $genomic_start = $block->start - $gap_lengths{5};
          }
			 else{#Just set to first start
			   $genomic_start = $block->start;
			 }
		   }
		  else{#Must be 5' Gap/overhang
			$gap_lengths{5} = $block->length;
		  }
		}
		elsif($block->isa('Bio::EnsEMBL::Mapper::Coordinate')){
		  $genomic_end = $block->end;
		}
		else{#Must be 3' Gap/overhang
		  $gap_lengths{3}  = $block->length;
		  $genomic_end    += $gap_lengths{3};
		}
		

		#So calculate gap coordinates
		#We want to skip and Gap objects
		#as these will have already been inserted into the cigarline by ExonerateProbe???
		#These are 5'/3' overhangs
		next if $block->isa('Bio::EnsEMBL::Mapper::Gap');

		if(@gaps){
		  push @gaps, ($block->start - $gaps[$#gaps]);
		}

		push @gaps, ($block->end + 1);
	  }

	  #remove last value as this is the end of the match and not a gap start
	  pop @gaps;

	  
	  #Insert intron gaps into cDNA cigarline
	  $gap_start  = shift @gaps;
	  $gap_length = shift @gaps;

	  #Set this to a hypothetical previous block end
	  #So the initial block_start calculation will be valid
	  $block_end  = ($genomic_start - 1);
	  my $last_gap       = 0;	
	  my $match_count    = 0;
	  my $mismatch_count = 0;
	  my $score          = 0;
	  my $q_length       = 0;


	  #Account for 5'/3' overhangs here by simple bp matching
	  #We need to account for Gap length as we may get 2m
	  #Which might actually be 1m1M, where the 1M is an overhang Gap
	  my ($gap_block);

	  foreach my $end(keys %gap_lengths){
		
		if($gap_lengths{$end}){
          #$gap_block here is actually mismatch & gap_block
          
		  if($end == 5){
			$gap_block = shift @stranded_cigar_line;
		  }
		  else{
			$gap_block = pop @stranded_cigar_line;
		  }
		  
          
		  @tmp = split//, $gap_block;
		  $align_type = pop @tmp;
		  ($align_length = $gap_block) =~ s/$align_type//;

          
          #Getting problem here with block_length is 1 but cigar is 2X
          #X > 1 here can span the end of the UTR and then an overhang!
          #dependant on how many mismatches are allowed

          #split this into two cigar blocks to represent the mismatch and the overhang
          #and replace the mismatched align block for processing below.
          

		  if( $align_type ne 'X' ) {
			throw("${end}' overhanging Gap has non-X/unexpected alignment type:\t$align_type");
		  }
          
       
          #if( $align_length > $gap_lengths{$end}){
			#throw("${end}' overhanging Gap has mismatch between cigar length".
            #      " ($align_length) and Gap->block_length".$gap_lengths{$end});
          #we have a mismatch cigar block which spans an alignment and a 5/3' gap/overhang
          #This was always handled below, but never got there because of this thrwo
          #}


		  
		  #Now match against genomic sequence???
		  #No way of doing this as we don't have the probe sequence!!!!
          #Can we get the query seq from ExonerateProbe?
		  $gap_block = $gap_lengths{$end}.'S';
		  $align_length -= $gap_lengths{$end};

      
		  if($align_length){
            #We have a mismatch cigar block which spans an alignment and a 5/3' gap/overhang
            #Replace mismatch align block

			if($end == 5 ){
              unshift @stranded_cigar_line, $align_length.$align_type;
			}
            else{ # is 3
              push @stranded_cigar_line, $align_length.$align_type;
			}
		  }

          #Add gap/overhang blocks
		  if($end == 5){
			unshift @stranded_cigar_line, $gap_block;
          }
          else{
			push @stranded_cigar_line, $gap_block;
		  }
		}
	  }

	  foreach my $block(@stranded_cigar_line){
		@tmp = split//, $block;
		$align_type = pop @tmp;
		($align_length = $block) =~ s/$align_type//;

		$q_length += $align_length;

		if($align_type eq '='){
		  $match_count += $align_length;
		  $score       += ($align_length * 5);
		}
		else{# 'X' or S  mismatch (was U)
		  $mismatch_count += $align_length;
		  $score          -= ($align_length * 4);
		}


		#We need to calculate the genomic end of the current block
		#given the previous block end or the genomic_start
		$block_start = ($block_end + 1);
		$block_end += $align_length;

		if($block_end >= $gap_start){

		  #Could have multiple deletions
		  while(($block_end >= $gap_start) && ! $last_gap){
			#Insert the match first
			$align_length = ($gap_start - $block_start);
			$cigar_line  .= $align_length.$align_type if $align_length;

			#Deletion wrt probe i.e. intron
			$cigar_line .= $gap_length.'D';

			#Now redefine start and end values
			#warn "block_start += $align_length + $gap_length";
			$block_start += $align_length + $gap_length;
			$block_end   += $gap_length;
		
			#Now grab the next gap
			if(@gaps){
			  $gap_start  = shift @gaps;
			  $gap_length =  shift @gaps;
			}
			else{
			  $last_gap = 1;
			}
		  }
		  
		  #We have reached the end of the gaps in this block
		  #so just redefine block here
		  $align_length = ($block_end - $block_start + 1);
		  $block = ($align_length) ? $align_length.$align_type : '';
		}

		$cigar_line .= $block;
	  }


	  #We could assign the start end directly
	  $feature->start($genomic_start);
	  $feature->end($genomic_end);
	  $feature->strand($transcript_strand);
	  $feature->cigar_string($cigar_line);
	  #warn "Final Feature $genomic_start $genomic_end $cigar_line" if $warn;


	  #Test if we have already seen this alignment
	  $gene       = $gene_adaptor->fetch_by_transcript_stable_id($seq_id);
	  #The only way of doing this is to test the genomic_start/end and the genomic cigarline with the gene_stable_id and the probe_id
	  $gene_sid = $gene->stable_id;
	  $gene_hit_key = "${gene_sid}:${probe_id}:${genomic_start}:${genomic_end}:${cigar_line}";


	  if(exists $gene_hits{$gene_hit_key}){
		$load_feature = 0;
	  }else{
		#No need to count hits here
		$gene_hits{$gene_hit_key} = undef;
	  }

	  #Now store the IDXref for this probe transcript hit
	  #This will mean we don't have recalculate this during the probe2transcript mapping step


	  #% ID over aligned region
	  #Which can be different if we have I|Ds
	  #As these will give different seq lengths
	  #But this is the ungapped alignment
	  #So ignore Target ID?
	  #cigar_line is reported wrt to strand of target
	  #we filter out all -ve hits by this point so don't need to account for here.
	  $query_perc = ($match_count/($match_count + $mismatch_count)) * 100;
	  $display_name = $self->get_display_name_by_stable_id($seq_id, 'transcript');

	  #$id_xref = Bio::EnsEMBL::IdentityXref->new
	  $xref = Bio::EnsEMBL::DBEntry->new
		(
		 #-XREF_IDENTITY => $query_perc,
		 #-TARGET_IDENTITY => 90.1,
		 #-SCORE => $score,
		 #-EVALUE => 12,
		 #-CIGAR_LINE => $cigar_line,
		 #-XREF_START => 1,#We are currently padding with mismatches to full length of query
		 #-XREF_END => $q_length,
		 #-ENSEMBL_START => $transcript_start,#target/hit_start
		 #-ENSEMBL_END => $transcript_end,#target/hit_end
		 #-ANALYSIS => $analysis,
		 -ANALYSIS => $feature->analysis,
		 -PRIMARY_ID => $seq_id,
		 -DISPLAY_ID => $display_name,
		 -DBNAME  => $edb_name,
		 -release => $schema_build,
		 -info_type => 'MISC',
		 -info_text => 'TRANSCRIPT',
		 -linkage_annotation => "ProbeTranscriptAlign $query_perc",#Add query_perc here when we have analysis
		 #-info_text => , #? What is this for? Is used in unique key so we get duplicated if null!!!
		 -version => $transcript_cache{$seq_id}->version, #version of transcript sid?

		);
	  #No strand here! Always +ve?!
	  #$dbe_adaptor->store($id_xref, $probe_id, 'Probe', 1);#Do we need ignore release flag here?


	  ###This cannot be done until we have ox.analysis_id in place v54?
	  #Yes we can, just use the linkage_annotation for now!


	  #No store this as a ProbeFeature DBEntry in line with how the individual probe xrefs
	  #are stored in probe2transcript
	  #we don't really need this extra translation and score info here
	  #and we are storing the cigar_line as part of the probe_feature
	  #so we can reconstitute the alignment if required

	  #This will mean we can separate the individual feature xrefs from the Probe/ProbeSet level full xref annotation
	  #Will mean some duplication for single probe arrays
	  #But will allow different sets of rules between ProbeAlign and probe2transcript
	  #i.e. we could relax the mapping rules but keep the xreffin rules more stringent

	  #We don't know what the ProbeFeature dbID is yet so we will have to pass this to store with the feature

	  #$dbe_adaptor->store($xref, $probe_id, 'Probe', 1);#Do we need ignore release flag here?
	}
	else{#No introns - reformat cigar line
	  $feature->cigar_string(join('', (split/:/, $feature->cigar_string)));
	}

	if($load_feature){

	  #Reset start ends for non-ref slices

	  if($slice->start != 1){

		#Need to reslice, but we don't want to affect
		#the slice in the cahse as this will screw up
		#further feature start/end tranforms

		my ($level, undef, $name) = split /:/, $slice->name;

		#warn "Resetting slice for ".$slice->name;
		my $start = $feature->start + $slice->start - 1;
		my $end   = $feature->end   + $slice->start - 1;
		$feature->start($start);
		$feature->end($end);
	
		warn "fetching slice $level,  $name, 1, $end"; 
	
		$slice = $slice_adaptor->fetch_by_region($level, $name, 1, $end);
	  }

	  #warn "Final start stop ".$feature->start.' '.$feature->end;

	  $feature->slice($slice);

	  # Recover the probe from the cache or DB
	  my $real_probe = $self->probes->{$probe_id};
	  
	  if(!$real_probe){

		$real_probe = $probe_adaptor->fetch_by_dbID($probe_id);
      
		if (!$real_probe){
		  throw "Trying to clean up features for persistence: cant find probe with dbID $probe_id in database!\n";
		}
		
		$self->probes->{$probe_id} = $real_probe;
	  }

	  $feature->probe($real_probe);
# 	  $feature->analysis($analysis);
	  push @features, [ $feature, $xref ];
	}
  }


  print "Finished set_probe_and_slice\n";

  return $self->features(\@features);
}


sub outdb {
  my ($self) = @_;

  my $outdb;
  my $dnadb;

  #Do we need to alter this???????????????????????????????????????????????????????????????????????????????????????????????????
  #Look at ImportArrays
  #There is method duplication here which we could move to a FuncgenDB.pm?
  #Do this when we move to hive


  if(! defined $self->{'efg_db'}){

	#This DNADB testing is a work around to avoid having to edit Config::ProbeAlign
	my $dnadb;

	#not defined as an empty env var will give the defined null string

	if($self->DNADB->{-dbname}){
	  $dnadb =  new Bio::EnsEMBL::DBSQL::DBAdaptor(%{ $self->DNADB });
	  $dnadb->dbc->disconnect_when_inactive(1);
	}

	  my $dba_hash = $self->OUTDB;
	  
	  my %dba_hash_copy = %$dba_hash;
	  
	  delete $dba_hash_copy{'-dnadb_name'};
	  delete $dba_hash_copy{'-dnadb_host'};
	  delete $dba_hash_copy{'-dnadb_port'};
	  delete $dba_hash_copy{'-dnadb_user'};
	
	$self->{'efg_db'} = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
	  (
	   %dba_hash_copy,
	   -dnadb => $dnadb,
	  );
	  
	  $self->{'efg_db'}->dbc->disconnect_when_inactive(1);
	  
   
	if(! $self->DNADB->{-dbname}){
	  print "WARNING: Using default DNADB ". $self->{'efg_db'}->dnadb->dbname."\n";
	}

	#We are now forcing use of OUTDB/DNADB config
	#else {
	#  #Historical, but sensible?
	#  $self->{'efg_db'} = $self->SUPER::db;
	#}
  }
  
  return $self->{'efg_db'};
}

sub query_file {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_query_file'} = $value;
  }

  if ( exists( $self->{'_query_file'} ) ) {
    return $self->{'_query_file'};
  } else {
    return undef;
  }
}


#############################################################

sub unmapped_objects{
  my ( $self, $value ) = @_;

  $self->{'_unmnapped_objects'} = $value if defined $value;
  
  return $self->{'_unmapped_objects'};
  
}

#############################################################


#############################################################

sub arrays {
  my ( $self, $value ) = @_;

  $self->{'_arrays'} = $value if defined $value;
  
  return $self->{'_arrays'};
  
}

#############################################################

sub probes {
  my ( $self, $value ) = @_;

  $self->{'_probes'} = $value if  defined $value;
  return $self->{'_probes'};
}

#############################################################

sub features {
  my ( $self, $value ) = @_;

  $self->{'_features'} = $value  if defined $value;
  return $self->{'_features'};#do we need to test and return undef?
}

sub analysis{
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    throw("Must pass RunnableDB:analysis a Bio::EnsEMBL::Analysis".
          "not a ".$analysis) unless($analysis->isa
                                     ('Bio::EnsEMBL::Analysis'));
    $self->{'analysis'} = $analysis;
  }
  return $self->{'analysis'};
}

sub output{
  my ($self, $output) = @_;
  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if($output){
    if(ref($output) ne 'ARRAY'){
      throw('Must pass RunnableDB:output an array ref not a '.$output);
    }
    push(@{$self->{'output'}}, @$output);
  }
  return $self->{'output'};
}

sub QUERYSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_QUERYSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_QUERYSEQS'} ) ) {
    return $self->{'_CONFIG_QUERYSEQS'};
  } else {
    return undef;
  }
}

sub QUERYTYPE {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_QUERYTYPE'} = $value;
  }

  if ( exists( $self->{'_CONFIG_QUERYTYPE'} ) ) {
    return $self->{'_CONFIG_QUERYTYPE'};
  } else {
    return undef;
  }
}

sub TARGETSEQS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_TARGETSEQS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_TARGETSEQS'} ) ) {
    return $self->{'_CONFIG_TARGETSEQS'};
  } else {
    return undef;
  }
}

sub MAX_MISMATCHES {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_MAX_MISMATCHES'} = $value;
  }

  if ( exists( $self->{'_CONFIG_MAX_MISMATCHES'} ) ) {
    return $self->{'_CONFIG_MAX_MISMATCHES'};
  } else {
    return undef;
  }
}

sub IIDREGEXP {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_IIDREGEXP'} = $value;
  }

  if ( exists( $self->{'_CONFIG_IIDREGEXP'} ) ) {
    return $self->{'_CONFIG_IIDREGEXP'};
  } else {
    return undef;
  }
}

sub OUTDB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OUTDB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OUTDB'} ) ) {
    return $self->{'_CONFIG_OUTDB'};
  } else {
    return undef;
  }
}

sub DNADB {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_DNADB'} = $value;
  }

  if ( exists( $self->{'_CONFIG_DNADB'} ) ) {
    return $self->{'_CONFIG_DNADB'};
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ( $self, $value ) = @_;

  if ( defined $value ) {
    $self->{'_CONFIG_OPTIONS'} = $value;
  }

  if ( exists( $self->{'_CONFIG_OPTIONS'} ) ) {
    return $self->{'_CONFIG_OPTIONS'};
  } else {
    return undef;
  }
}

sub HIT_SATURATION_LEVEL {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{'_CONFIG_HIT_SATURATION_LEVEL'} = $val;
  }

  if (exists $self->{'_CONFIG_HIT_SATURATION_LEVEL'}) {
    return $self->{'_CONFIG_HIT_SATURATION_LEVEL'};
  } else {
    # default to 100
	# This is not used for ProbeTranscriptAlign
    return 100;
  }
}

###############################################
###     end of config
###############################################

1;
