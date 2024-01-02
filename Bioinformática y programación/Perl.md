
## Installing modules in local

[cpanminus](http://search.cpan.org/~miyagawa/App-cpanminus-1.4008/lib/App/cpanminus.pm) is quickly becoming the choice interface for CPAN. It supports installing packages in to the user's home directory.

```perl
curl -L http://cpanmin.us | perl - App::cpanminus
```

```perl
curl -L http://cpanmin.us | perl - Lingua::Romana::Perligata
```

```perl
export PERL5LIB=$HOME/perl5/lib/perl5:$PERL5LIB
```