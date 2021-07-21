# dunes

A command line tool for introducing random mutations into DNA sequences.

## Download or Compilation

- Download:
[dunes.jar](https://github.com/Sergey-Knyazev/dunes/releases/download/0.1/dunes.jar)

- Compilation:
```bash
mvn clean install
```

## Running

```bash
java -jar dunes.jar
```

## Help

```bash
java -jar dunes.jar -h

Usage: dunes [-hV] -i=FILE [-m=m] [-n=n] [-o=<outputFolder>] [-y=y]
  -h, --help               Show this help message and exit.
  -i, --inFile=FILE        input fasta file with viral sequences to be mutated
  -m, --mutation-rate=m    mutation rate in substitutions per nucleotide per year
                             (s/n/y) (default: 4.1e-3)
  -n, --mutants-number=n   number of mutants for a strain (default: 1)
  -o, --outFile=<outputFolder>
                           fasta file with mutated sequences (default: "<inFile
                             base>"Mut.fasta")
  -V, --version            Print version information and exit.
  -y, --years=y            years of evolution (default: 1)
```
