
#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor
#
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor - A database aware adaptor responsible for
the creation of Slices in the context of eFG objects.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
 

  $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);

  $slice_adaptor = $db->get_SliceAdaptor();

  # get a slice on the entire chromosome X
  $gene_regulation_slice = $slice_adaptor->fetch_by_gene($gene);


=head1 DESCRIPTION

This module is simple wrapper class for the core SliceAdaptor, extending new
methods to generate Slices for eFG features associated with a given gene or transcript.

=head1 CONTACT

This module is part of the Ensembl project http://www.ensembl.org

For more information email <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut


package Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;


use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning stack_trace_dump);


@ISA = ('Bio::EnsEMBL::DBSQL::SliceAdaptor');


=head2 fetch_by_Gene_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : my $gene_reg_slice = $efg_slice_adaptor->fetch_by_gene($gene);
  Description: Fetches a slice by finding the associated regulatory
               elements associated with a given gene(and it's transcripts and 
               translations), and extending the gene feature slice to 
               encapsulate the associated regulatory elements.
  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if incorrent arg provided
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_Gene_FeatureSets {
  my ($self, $gene, $fsets) = @_;

  if(! ( ref($gene) && $gene->isa('Bio::EnsEMBL::Gene'))){
    throw("You must pass a valid Bio::EnsEMBL::Gene object");
  }
  
  #Now we want to get all the eFG feature associated with this gene and build the slice
  #resitrict to RegFeats for now.
  my $cs_name      = $gene->slice->coord_system_name;
  my $gene_chr     = $gene->seq_region_name;
  my $start        = $gene->seq_region_start;
  my $end          = $gene->seq_region_end;
  

  ($start, $end) = $self->_set_bounds_by_xref_FeatureSets($gene_chr, $start, $end, $gene, $fsets);

  #Now need to do this for all transcripts and proteins too?
  
  foreach my $trans(@{$gene->get_all_Transcripts}){
	($start, $end) = $self->_set_bounds_by_xref_FeatureSets($gene_chr, $start, $end, $trans, $fsets);
	($start, $end) = $self->_set_bounds_by_xref_FeatureSets($gene_chr, $start, $end, $trans->translation, $fsets);
  }

  return $self->fetch_by_region($cs_name, $gene_chr, $start, $end, $gene->strand);
}

=head2 fetch_by_Transcript_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : my $trans_reg_slice = $efg_slice_adaptor->fetch_by_transcript($transcript);
  Description: Fetches a slice by finding the associated regulatory
               elements associated with a given transcript(and it's translation), and 
               extending the gene feature slice to encapsulate the associated regulatory 
               elements.
  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if incorrent arg provided
  Caller     : self
  Status     : at risk

=cut

sub fetch_by_Transcript_FeatureSets{
  my ($self, $transcript, $fsets) = @_;


  if(! ( ref($transcript) && $transcript->isa('Bio::EnsEMBL::Transcript'))){
    throw("You must pass a valid Bio::EnsEMBL::Transcript object");
  }
  
  #Now we want to get all the eFG feature associated with this gene and build the slice
  #resitrict to RegFeats for now.
  my $cs_name   = $transcript->slice->coord_system_name;
  my $trans_chr = $transcript->seq_region_name;
  my $start     = $transcript->seq_region_start;
  my $end       = $transcript->seq_region_end;
  

  ($start, $end) = $self->_set_bounds_by_xref_FeatureSets($trans_chr, $start, $end, $transcript, $fsets);

  #Now need to do this for the protein too
 
  ($start, $end) = $self->_set_bounds_by_xref_FeatureSets($trans_chr, $start, $end, $transcript->translation, $fsets);



  return $self->fetch_by_region($cs_name, $trans_chr, $start, $end, $transcript->strand);
}

#Do we need fetch_by_translation?

=head2 _set_bounds_by_xref_FeatureSets

  Arg [1]    : string - seq_region_name i.e. chromosome name.
  Arg [2]    : string - seq_region_start of current slice bound
  Arg [3]    : string - seq_region_end of current slice bound.
  Arg [4]    : Bio::EnsEMBL::Gene|Transcript|Translation
  Arg [5]    : arrayref - Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : ($start, $end) = $self->_set_bounds_by_regulatory_feature_xref
                              ($trans_chr, $start, $end, $transcript, $fsets);
  Description: Internal method to set an xref Slice bounds given a list of
               FeatureSets.
  Returntype : List - ($start, $end);
  Exceptions : throw if incorrent args provided
  Caller     : self
  Status     : at risk

=cut

sub _set_bounds_by_xref_FeatureSets{
  my ($self, $chr, $start, $end, $obj, $fsets) = @_;

  my ($extdb_name, $efg_feature);
  my $dbe_adaptor = $self->db->get_DBEntryAdaptor;


  #do we need to test for start/end/chr here?
  #Is implicit if we test for obj and fsets, but 
  #does not check if defined.

  #Set ext_dbname and validate obj
  #Do we need a central store for ensembl db names?

  if($obj->isa('Bio::EnsEMBL::Gene')){
	$extdb_name = 'ensembl_core_Gene';
  }
  elsif($obj->isa('Bio::EnsEMBL::Transcript')){
	$extdb_name = 'ensembl_core_Transcript';
  }
  elsif($obj->isa('Bio::EnsEMBL::Translation')){
	$extdb_name = 'ensembl_core_Translation';
  }
  else{
	throw('Currently only handles Ensembl Gene, Transcript and Translation xrefs');
  }


  #warn "Setting bounds by $obj ".$obj->stable_id;

  #Set which eFG features we want to look at.

  if(! (defined $fsets || ( ref($fsets) eq 'ARRAY' && scalar(@$fsets) == 0)) ){
	throw('Must define an array of Bio::EnsEMBL::FeatureSets to extend xref Slice bound by');
  }

  my %feature_set_types;

  foreach my $fset(@$fsets){
	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
	
	$feature_set_types{$fset->type} ||= [];
	push @{$feature_set_types{$fset->type}}, $fset;
  }


  #Get xrefs for each eFG feature set type
  foreach my $fset_type(keys %feature_set_types){

	my $xref_method    = 'list_'.$fset_type.'_feature_ids_by_extid';
	#e.g. list_regulatory_feature_ids_by_extid

	my $adaptor_method = 'get_'.ucfirst($fset_type).'FeatureAdaptor';

	my $adaptor = $self->db->$adaptor_method;

	my %feature_set_ids;
	map $feature_set_ids{$_->dbID} = 1, @{$feature_set_types{$fset_type}};

	my $cnt = 0;

	foreach my $dbID($dbe_adaptor->$xref_method($obj->stable_id, $extdb_name)){
	  $efg_feature = $adaptor->fetch_by_dbID($dbID);
	
	  #Skip if it is not in one of our FeatureSets
	  next if ! exists $feature_set_ids{$efg_feature->feature_set->dbID};
	  next if $efg_feature->seq_region_name ne $chr;
	  

	  #warn "found xref ".$efg_feature->display_label.' with start '.$efg_feature->seq_region_start;
	  $cnt ++;
  	  
	  $start = $efg_feature->seq_region_start if $efg_feature->seq_region_start < $start;
	  $end   = $efg_feature->seq_region_end   if $efg_feature->seq_region_end   > $end;
	}
	
	#	warn "Found $cnt $fset_type xrefs";

  }



  return ($start, $end);
}


1;
