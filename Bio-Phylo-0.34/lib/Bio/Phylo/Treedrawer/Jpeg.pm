package Bio::Phylo::Treedrawer::Jpeg;
use strict;
use Bio::Phylo::Util::Exceptions 'throw';
use vars qw(@ISA);

eval { require Bio::Phylo::Treedrawer::Png };
if ( $@ ) {
    throw 'ExtensionError' => "Error loading the Bio::Phylo::Treedrawer::Png extension: $@";
}
@ISA=qw(Bio::Phylo::Treedrawer::Png);

=head1 NAME

Bio::Phylo::Treedrawer::Jpeg - Graphics format writer used by treedrawer, no
serviceable parts inside

=head1 DESCRIPTION

This module creates a jpeg file from a Bio::Phylo::Forest::DrawTree
object. It is called by the L<Bio::Phylo::Treedrawer> object, so look there to
learn how to create tree drawings.


=begin comment

# only need to override finish to write to a different format

=end comment

=cut

sub _finish {
    my $self = shift;
    my $jpg;
    eval { $jpg = $self->_api->jpeg };
    if ( not $@ ) {
        return $jpg;
    }
    else {
        throw 'ExtensionError' => "Can't create JPEG, libgd probably not compiled with it"
    }
}

=head1 SEE ALSO

=over

=item L<Bio::Phylo::Treedrawer>

The jpeg treedrawer is called by the L<Bio::Phylo::Treedrawer> object. Look there
to learn how to create tree drawings.

=item L<Bio::Phylo::Manual>

Also see the manual: L<Bio::Phylo::Manual> and L<http://rutgervos.blogspot.com>.

=back

=head1 REVISION

 $Id: Jpeg.pm 1401 2010-10-22 15:02:22Z rvos $

=cut

1;





