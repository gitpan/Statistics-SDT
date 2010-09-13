use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::SDT',
    license             => 'perl',
    dist_author         => 'Roderick Garton <rgarton@cpan.org>',
    dist_version_from   => 'lib/Statistics/SDT.pm',
    build_requires => {
        'Test::More' => 0,
        'Math::Cephes' => 0.44,
    },
    add_to_cleanup      => [ 'Statistics-SDT-0.40-*' ],
);

$builder->create_build_script();