use strict;
use warnings;
use Test::More tests => 20;
BEGIN { use_ok('Statistics::SDT') };

my $sdt = Statistics::SDT->new(
	correction => 1,
  	precision_s => 2,
);
isa_ok($sdt, 'Statistics::SDT');

eval {
    $sdt->init(
      hits => 50,
  	signal_trials => 50,
	false_alarms => 17,
  	noise_trials => 25,
    );
};
ok(!$@);


my %refvals = (
	d => 1.86,
	Ad => .91,
	Ap => .82,
    beta => .07,
    logbeta => -2.6,
    cdecision => -1.4,
    griers => -.91,
    criterion => -.47,
    hr => .99,
    far => .68,
);

my $v;

foreach (qw/hr far/) {
    ok(defined $sdt->{$_} );
    ok(equal($sdt->{$_}, $refvals{$_}));
}

foreach (qw/d Ad Ap/) {
	$v = $sdt->sensitivity($_);
    ok( equal($v, $refvals{$_}), "Sensitivity $_ : $v = $refvals{$_}");
}

foreach (qw/beta logbeta cdecision griers/) {
	$v = $sdt->bias($_);
    ok( equal($v, $refvals{$_}), "Bias $_ $v = $refvals{$_}");
}

$v = $sdt->criterion();
ok( equal($v, $refvals{'criterion'}), "Criterion k $v = $refvals{'criterion'}");

$v = $sdt->dc2hr();
ok( equal($v, $refvals{'hr'}), "Hr $v = $refvals{'hr'}");

$v = $sdt->dc2far();
ok( equal($v, $refvals{'far'}), "FAr $v = $refvals{'far'}");

$v = $sdt->dc2logbeta();
ok( equal($v, $refvals{'logbeta'}), "Logbeta $v = $refvals{'logbeta'}");

$v = $sdt->sensitivity('f' => {hr => .866, states => 3, correction => 0, method => 'alexander'});
ok( equal($v, 2), "Sensitivity fc $v = 2");

$v = $sdt->sensitivity('f' => {hr => .866, states => 3, correction => 0, method => 'smith'});
ok( equal($v, 2.05), "Sensitivity fc $v = 2.05");

sub equal {
    return $_[0] == $_[1] ? 1 : 0;
}
