#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor - A database adaptor for fetching and
storing Funcgen feature sets.

=head1 SYNOPSIS

my $fs_adaptor = $db->get_FeatureSetAdaptor();


my @displayable_fsets = @{$fs_adaptor->fetch_all_displayable()};

=head1 DESCRIPTION

The FeatureSetAdaptor is a database adaptor for storing and retrieving
Funcgen feature set.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor; #DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor);

=head2 fetch_all_by_feature_class

  Arg [1]    : String - feature class i.e. 'annotated', 'regulatory', 'segmentation', 'mirna' or 'external'
  Arg [2]    : String (optional) - status e.g. 'DISPLAYABLE'
  Arg [2]    : Bio::EnsEMBL::Funcgen::Epigenome (optional) or a HASH parameters
               containing contraint config e.g.

                   $feature_set_adaptor->fetch_all_by_feature_class
                                           ('annotated',
                                             {'constraints' =>
                                               {
                                               epigenomes     => [$epigenome], #Bio::EnsEMBL::Funcgen::Epigenome
                                               projects       => ['ENCODE'],
                                               evidence_types => ['Hists & Pols'],
                                               feature_types  => [$ftype], #Bio::EnsEMBL::Funcgen::FeatureType
                                               }
                                             });

  Example    : my @fsets = @{$fs_adaptopr->fetch_all_by_feature_class('annotated')};
  Description: Retrieves FeatureSet objects from the database based on feature_set type.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if type not defined
  Caller     : General
  Status     : At Risk - Move status to params hash

=cut

sub fetch_all_by_feature_class {
  my ($self, $type, $status, $params) = @_;

  throw('Must provide a feature_set type') if(! defined $type);
  my $sql = "fs.type = '".$type."'";

  if (defined $params) {        #Some redundancy over $epigenome arg and $params epigenome

    if ( ref($params) eq 'Bio::EnsEMBL::Funcgen::Epigenome') {
      $params = {constraints => {epigenomes => [$params]}};
    } elsif (ref($params) ne 'HASH') {
      throw('Argument must be a Bio::EnsEMBL::Funcgen::Epigenome or a params HASH');
    }
  }


  if ($status) {
    $params->{constraints}{states} = [ $status ];
  }

  #Deal with params constraints
  my $constraint = $self->compose_constraint_query($params);
  $sql .=  " AND $constraint " if $constraint;


  #Get result and reset true tables
  my $result = (defined $sql) ? $self->generic_fetch($sql) : [];
  $self->reset_true_tables;

  return $result;
}


=head2 fetch_all_displayable_by_type

  Arg [1]    : String - Type of feature set i.e. 'annotated', 'regulatory' or 'supporting'
  Arg [2]    : Bio::EnsEMBL::Funcgen::Epigenome (optional) or parameters HASH
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_type('annotated');
  Description: Wrapper method for fetch_all_by_type
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_displayable_by_type {
    my ($self, $type, $epigenome_or_params) = @_;

	#Move status to config hash
	$self->fetch_all_by_feature_class($type, 'DISPLAYABLE', $epigenome_or_params);

}


=head2 fetch_by_name

  Arg [1]    : string - name of FeatureSet
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = @{$fset_adaptor->fetch_by_name('feature_set-1')};
  Description: Fetch all FeatureSets wit a given name
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if no name passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name {
  my ($self, $name, $status) = @_;

  throw("Must provide a name argument") if (! defined $name);

  my $sql = "fs.name='".$name."'";

  if($status){
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return (defined $sql) ? $self->generic_fetch($sql)->[0] : [];

}


=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  return ([ 'feature_set', 'fs' ]);
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
	my $self = shift;

	return qw( fs.feature_set_id fs.feature_type_id
			   fs.analysis_id fs.epigenome_id
			   fs.name fs.type
			   fs.description fs.display_label
			   fs.experiment_id);
}




=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@fsets, $fset, $analysis, %analysis_hash, $feature_type, $epigenome, $name, $type, $display_label, $desc);
	my ($feature_set_id, $ftype_id, $analysis_id, $epigenome_id, $exp_id, %ftype_hash, %epigenome_hash);

	my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
	my $anal_adaptor = $self->db->get_AnalysisAdaptor();
	my $epigenome_adaptor = $self->db->get_EpigenomeAdaptor();
	$epigenome_hash{'NULL'} = undef;

	$sth->bind_columns(\$feature_set_id, \$ftype_id, \$analysis_id, \$epigenome_id,
                       \$name, \$type, \$desc, \$display_label, \$exp_id);

	while ( $sth->fetch()) {
		$epigenome_id ||= 'NULL';

		# Get the analysis object
		$analysis_hash{$analysis_id} = $anal_adaptor->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});

		# Get the feature type object
		$ftype_hash{$ftype_id} = $ft_adaptor->fetch_by_dbID($ftype_id) if(! exists $ftype_hash{$ftype_id});

		# Get the epigenome object
        $epigenome_hash{$epigenome_id}
            = $epigenome_adaptor->fetch_by_dbID($epigenome_id)
            if ( !exists $epigenome_hash{$epigenome_id} );


		#Use new_fast here and strip the prefixed -'s
		$fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
		  (
		   -dbID          => $feature_set_id,
		   -adaptor       => $self,
		   -feature_type  => $ftype_hash{$ftype_id},
		   -analysis      => $analysis_hash{$analysis_id},
		   -epigenome     => $epigenome_hash{$epigenome_id},
		   -name          => $name,
		   -feature_class => $type,
		   -display_label => $display_label,
		   -description   => $desc,
		   -experiment_id => $exp_id,
		  );

		push @fsets, $fset;

	}

	return \@fsets;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Example    : $oaa->store($fset1, $fset2, $fset3);
  Description: Stores FeatureSet objects in the database.
  Returntype : Listref of stored FeatureSet objects
  Exceptions : Throws if FeatureSet does not have a stored FeatureType
               Throws if invalid FeatureSet passed
               Throws if not FeatureSets passed
               Warns if external_db_name not defined is type is external
               Throws if external_db is not present in the db
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @fsets = @_;

	throw('Must supply a list of FeatureSets to store') if(scalar(@fsets) == 0);

	my $sth = $self->prepare
	  (
	   "INSERT INTO feature_set
        (feature_type_id, analysis_id, epigenome_id, name, type, description, display_label, experiment_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
	  );

	my $db = $self->db;
	my ($sql, $edb_id, %edb_hash);

    foreach my $fset (@fsets) {
		throw('Can only store FeatureSet objects, skipping $fset')	if ( ! $fset->isa('Bio::EnsEMBL::Funcgen::FeatureSet'));


		if (! $fset->is_stored($db) ) {

		  # Check FeatureType and Analysis
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fset->feature_type);
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $fset->analysis);


		  # Check optional Epigenome and Experiment
		  my $epigenome_id;
		  my $epigenome = $fset->epigenome;

		  if($epigenome){
			$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Epigenome', $epigenome);
			$epigenome_id = $epigenome->dbID;
		  }

		  my $exp_id;
		  my $exp = $fset->experiment;

		  if($exp){
        $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exp);
        $exp_id = $exp->dbID;
		  }


		  $sth->bind_param(1, $fset->feature_type->dbID, SQL_INTEGER);
		  $sth->bind_param(2, $fset->analysis->dbID,     SQL_INTEGER);
		  $sth->bind_param(3, $epigenome_id,             SQL_INTEGER);
		  $sth->bind_param(4, $fset->name,               SQL_VARCHAR);
		  $sth->bind_param(5, $fset->feature_class,      SQL_VARCHAR);
		  $sth->bind_param(6, $fset->description,        SQL_VARCHAR);
		  $sth->bind_param(7, $fset->display_label,      SQL_VARCHAR);
		  $sth->bind_param(8, $exp_id,                   SQL_INTEGER);

		  $sth->execute();
		  $fset->dbID($self->last_insert_id);
		  $fset->adaptor($self);
		}
		else{
			warn('FeatureSet '.$fset->name.'is already stored, updating status entries');
			$self->store_states($fset);
		}
	}
	return \@fsets;
}

#focus is old terminology
#tend to use 'core' now
#although here we use focus to determine a set which is actually used in the build
#where as core is at the feature type level
#i.e. you can have a feature set with a core feature type
#which is not in the build and is therefore not a focus set.

=head2 fetch_focus_set_config_by_FeatureSet

  Args       : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : $self->{'focus_set'} = $self->adaptor->fetch_focus_set_config_by_FeatureSet($self);
  Description: Caches and returns focus set config for a given FeatureSet
  Returntype : Boolean
  Exceptions : Warns if meta entry not present
  Caller     : Bio::EnsEMBL::Funcgen::FeatureSet::is_focus_set
  Status     : At Risk

=cut

sub fetch_focus_set_config_by_FeatureSet{
  my ($self, $fset) = @_;

  $self->{focus_set_config} ||= {};

  if (! defined $self->{focus_set_config}->{$fset->dbID}) {

    #Is is an attribute set?
    if ($self->fetch_attribute_set_config_by_FeatureSet($fset)) {

      #Need to define these as RegBuild config
      if ($fset->feature_type->is_core_evidence) {
        $self->{focus_set_config}->{$fset->dbID} = 1;
      }
    }
  }

  return $self->{focus_set_config}->{$fset->dbID};
}


=head2 fetch_attribute_set_config_by_FeatureSet

  Args       : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : $self->{'attribute_set'} = $self->adaptor->fetch_attribute_set_config_by_FeatureSet($self);
  Description: Caches and returns attribute set config for a given FeatureSet
  Returntype : Boolean
  Exceptions : Warns if meta entry not present
  Caller     : Bio::EnsEMBL::Funcgen::FeatureSet::is_attribute_set
  Status     : At Risk

=cut

sub fetch_attribute_set_config_by_FeatureSet{
   deprecate(
    "Bio::EnsEMBL::Funcgen::FeatureSet::fetch_attribute_set_config_by_FeatureSet has been deprecated and will be removed in Ensembl release 89.");
}


sub fetch_feature_set_filter_counts{
  my $self = shift;

   my $sql = 'SELECT count(*), eg.name, eg.description, eg.is_project, ft.class, epi.name, epi.description '.
    'FROM experimental_group eg, experiment e, feature_set fs, feature_type ft, epigenome epi '.
        'WHERE fs.experiment_id=e.experiment_id '.
          'AND e.experimental_group_id=eg.experimental_group_id '.
          'AND fs.feature_type_id=ft.feature_type_id AND fs.epigenome_id=epi.epigenome_id '.
            'AND fs.type="annotated" '.
                'GROUP BY eg.name, eg.is_project, ft.class, epi.name';

  #warn $sql;
  #Need to write HC around this as we sometimes get less than expect.


  my @rows = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
  my $ft_a = $self->db->get_FeatureTypeAdaptor;
  my $ftype_info = $ft_a->get_regulatory_evidence_info;

  my %filter_info =
    (
     #Project=> {},
     #'Cell/Tissue' => {},
     All =>
     {
      All =>{ count       => 0,
              description => 'All experiments',
            }
     }

    );

  foreach my $row(@rows){

    my ($count, $project, $proj_desc, $is_proj, $ft_class, $epigenome_name, $epigenome_desc) = @$row;

    #All counts
    $filter_info{All}{All}{count} += $count;

    #Project counts
    if($is_proj){

      if(! exists $filter_info{Project}{$project}){
        $filter_info{Project}{$project} =
          { count       => 0,
            description => $proj_desc,
          };
      }

      $filter_info{Project}{$project}{count} += $count;
    }

    #Cell/Tissue counts
    if(! exists $filter_info{'Cell/Tissue'}{$epigenome_name}){
      $filter_info{'Cell/Tissue'}{$epigenome_name} =
        { count       => 0,
          description => $epigenome_desc,
        };
    }
    $filter_info{'Cell/Tissue'}{$epigenome_name}{count} += $count;

    #Evidence class counts
    #Do we want to split this into ft.class
    #i.e. split 'DNase1 & TFBS'
    my $evidence_type  = $ft_a->get_regulatory_evidence_type($ft_class);
    my $ft_class_label =  $ftype_info->{$evidence_type}{label};

    if(! exists $filter_info{'Evidence type'}{$ft_class_label}){
      $filter_info{'Evidence type'}{$ft_class_label} =
        { count       => 0,
          description => $ftype_info->{$evidence_type}{long_name},
        };
    }
    $filter_info{'Evidence type'}{$ft_class_label}{count} += $count;
  }

  return \%filter_info;

  #Do we need to add an 'in_build' filter /data field?

}




#Need to bind param these as they come from URL parameters and are not tested

#Could move a lot of these to the BaseAdaptor if we have a valid_constraints config available
#and we use generic main_table approach for building constraint

#and reuse between adaptors if we use the _tables method to get the table syn
#This may mean contraints can be specified for classes which do not contain
#the relevant fields.
#Allow this flexiblity or validate fields/constraint?
#Or implicit by location of contraint config, i.e. put it in the relevant
#parent adaptors


#All these _constrain methods must return a valid constraint string, and a hashref of any other constraint config

sub _constrain_projects{
  return $_[0]->_constrain_experimental_groups($_[1]);
}


sub _constrain_experimental_groups{
  my ($self, $egs, $projects_only) = @_;

  #enable query extension
  my $constraint_conf = {tables => [['experiment', 'e']]};


  if ( (ref($egs) ne 'ARRAY') ||
       scalar(@$egs) == 0 ) {
    throw('Must pass an arrayref of project names');
  }
  my @eg_ids;

  foreach my $eg (@$egs) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ExperimentalGroup', $eg);

    if ( $projects_only && (! $eg->is_project) ) {
      throw("You have passed an ExperimentalGroup which is not a project:\t".$eg->name);
    }

    push @eg_ids, $eg->dbID;
  }

  #Don't need to bind param this as we validate above
  my $constraint = ' fs.experiment_id=e.experiment_id AND '.
    'e.experimental_group_id IN ('.join(', ', @eg_ids).')';

  return ($constraint, $constraint_conf);
}




sub _constrain_evidence_types {
  my ($self, $etypes) = @_;

	my $constraint_conf = {tables => [['feature_type', 'ft']]};

  my %in_values =
		(
		 'DNase1 & TFBS' => ['Open Chromatin',
                         'Transcription Factor',
                         'Transcription Factor Complex'],

		 'Hists & Pols'  => ['Histone',
                         'Polymerase'],
		);

  my @ft_classes;

  if ( (ref($etypes) ne 'ARRAY') ||
       scalar(@$etypes) == 0 ) {
		throw('Must pass an arrayref of evidence types');
  }

  foreach my $etype (@$etypes) {

		if (! exists $in_values{$etype}) {
		  throw("You have passed an invalid evidence type filter argument($etype)\n".
            "Please use one of the following:\t".join(' ,', keys(%in_values)));
		}
		push @ft_classes, @{$in_values{$etype}};
  }

  #Don't need to bind param this as we validate above
  my $constraint = ' fs.feature_type_id=ft.feature_type_id AND ft.class IN ("'.
		join('", "', @ft_classes).'")';

  return ($constraint, $constraint_conf);
}


sub fetch_all_by_type { # deprecated in v74, but still used by web
  my $self = shift;
  my $type = shift;
  my $status = shift;

  deprecate('Please use fetch_all_by_feature_class');

  return $self->fetch_all_by_feature_class($type, $status);
}

1;



