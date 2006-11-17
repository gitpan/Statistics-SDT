use 5.008008;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME              => 'Statistics::SDT',
    VERSION_FROM      => 'lib/Statistics/SDT.pm', # finds $VERSION
    PREREQ_PM         => { Math::Cephes => 0.43, Carp => 0 }, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Statistics/SDT.pm', # retrieve abstract from module
       AUTHOR         => 'Roderick Garton <rgarton@utas_._edu_._au>') : ()),
);
