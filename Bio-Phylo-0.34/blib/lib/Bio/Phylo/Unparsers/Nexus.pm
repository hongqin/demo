# $Id: Nexus.pm 1532 2010-11-26 17:20:35Z rvos $
package Bio::Phylo::Unparsers::Nexus;
use strict;
use Bio::Phylo::Unparsers::Abstract;
use Bio::Phylo::Util::CONSTANT qw(:objecttypes);
use Bio::Phylo::Util::Exceptions 'throw';
use vars qw(@ISA);
@ISA=qw(Bio::Phylo::Unparsers::Abstract);

=head1 NAME

Bio::Phylo::Unparsers::Nexus - Serializer used by Bio::Phylo::IO, no serviceable parts inside

=head1 DESCRIPTION

This module turns a L<Bio::Phylo::Matrices::Matrix> object into a nexus
formatted matrix. It is called by the L<Bio::Phylo::IO> facade, don't call it
directly. You can pass the following additional arguments to the unparse call:
	
	# an array reference of matrix, forest and taxa objects:
	-phylo => [ $block1, $block2 ]
	
	# the arguments that can be passed for matrix objects, 
	# refer to Bio::Phylo::Matrices::Matrix::to_nexus:
	-matrix_args => {}

	# the arguments that can be passed for forest objects, 
	# refer to Bio::Phylo::Forest::to_nexus:
	-forest_args => {}

	# the arguments that can be passed for taxa objects, 
	# refer to Bio::Phylo::Taxa::to_nexus:
	-taxa_args => {}	
	
	OR:
	
	# for backward compatibility:
	-phylo => $matrix	

=begin comment

 Type    : Wrapper
 Title   : _to_string($matrix)
 Usage   : $nexus->_to_string($matrix);
 Function: Stringifies a matrix object into
           a nexus formatted table.
 Alias   :
 Returns : SCALAR
 Args    : Bio::Phylo::Matrices::Matrix;

=end comment

=cut

sub _to_string {
    my $self   = shift;
    my $blocks = $self->{'PHYLO'};
    my $nexus  = "#NEXUS\n";
    my $type;
    eval { $type = $blocks->_type  };

    # array?
    if ( $@ ) {
    	for my $block ( @$blocks ) {
    		eval { $type = $block->_type };
    		my %args;
    		if ( $type == _FOREST_ ) {
    			if ( exists $self->{'FOREST_ARGS'} ) {
    				%args = %{ $self->{'FOREST_ARGS'} };
    			}
    		}
    		elsif ( $type == _TAXA_ ) {
    			if ( exists $self->{'TAXA_ARGS'} ) {
    				%args = %{ $self->{'TAXA_ARGS'} };
    			}    			
    		}
    		elsif ( $type == _MATRIX_ ) {
    			if ( exists $self->{'MATRIX_ARGS'} ) {
    				%args = %{ $self->{'MATRIX_ARGS'} };
    			}     			
    		}
    		elsif ( $@ ) {
		    	throw 'ObjectMismatch' => "Can't unparse this object: $blocks";
    		}
    		$nexus .= $block->to_nexus(%args);
    	}
    }
    
    # matrix?
    elsif ( defined $type and $type == _MATRIX_ ) {
    	$nexus .= $blocks->to_nexus;
    }
    
    # project?
    elsif ( defined $type and $type == _PROJECT_ ) {
    	$nexus = $blocks->to_nexus;
    }    
    
    # wrong!
    else {
    	throw 'ObjectMismatch' => "Can't unparse this object: $blocks";
    }
    
    return $nexus;

}

# podinherit_insert_token

=head1 SEE ALSO

=over

=item L<Bio::Phylo::IO>

The nexus serializer is called by the L<Bio::Phylo::IO> object.

=item L<Bio::Phylo::Manual>

Also see the manual: L<Bio::Phylo::Manual> and L<http://rutgervos.blogspot.com>.

=back

=head1 REVISION

 $Id: Nexus.pm 1532 2010-11-26 17:20:35Z rvos $

=cut

1;
