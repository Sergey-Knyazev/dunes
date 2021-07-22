import java.io.*;
import java.util.*;
import java.lang.Math;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import picocli.CommandLine;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet; 
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;

@CommandLine.Command(name = "dunes", mixinStandardHelpOptions = true, version = "0.0")
public class Main implements Runnable{
    @CommandLine.Option(names={"-i", "--inFile"}, description="input fasta file with viral sequences to be mutated",
            paramLabel = "FILE", required=true)
    private File inputFile;
    @CommandLine.Option(names={"-o", "--outFile"},
            description="fasta file with mutated sequences (default: \"<inFile base>\"Mut.fasta\")",
            paramLabel = "<outputFolder>")
    private File outputFile;
    @CommandLine.Option(names={"-m", "--mutation-rate"},
            description="mutation rate in substitutions per nucleotide per year (s/n/y) (default: 4.1e-3)",
            paramLabel = "m")
    private double rate = 4.1e-3;
    @CommandLine.Option(names={"-y", "--years"},
            description="years of evolution (default: 1)",
            paramLabel = "y")
    private double years = 1;
    @CommandLine.Option(names={"-n", "--mutants-number"},
            description="number of mutants for a strain (default: 1)",
            paramLabel = "n")
    private double n = 1;
    private double nuc_mut_prob;
    private static List<Character> nucs = new ArrayList<Character>() {{add('A'); add('C'); add('T'); add('G');}};
    private static Random rand = new Random();

    public void run() {
        try {
            nuc_mut_prob = Math.pow(1+rate,years) - 1;
            AmbiguityDNACompoundSet ambiguityDNACompoundSet = AmbiguityDNACompoundSet.getDNACompoundSet();
            DNASequenceCreator ambigDNASequenceCreator = new DNASequenceCreator(ambiguityDNACompoundSet);
            GenericFastaHeaderParser fastaHeaderParser = new GenericFastaHeaderParser();
            FastaReader fastaReader = new FastaReader(inputFile, fastaHeaderParser, ambigDNASequenceCreator);
            LinkedHashMap<String, DNASequence> seqs = fastaReader.process(); 
            Set<String> seq_names = seqs.keySet();
            LinkedHashSet<DNASequence> out_seqs = new LinkedHashSet<>();
            if(outputFile == null) {
                String[] tokens = inputFile.toString().split("\\.(?=[^\\.]+$)");
                outputFile = new File(String.format("%sMut.%s",tokens[0],tokens[1]));
            }
            for (String seq_name : seq_names) {
                 out_seqs.addAll(mutate_seqs(seqs.get(seq_name), seq_name));
            }
            FastaWriterHelper.writeNucleotideSequence(outputFile, out_seqs);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
    }

    private LinkedHashSet<DNASequence> mutate_seqs(DNASequence seq, String seq_name)
            throws CompoundNotFoundException {
        LinkedHashSet<DNASequence> out_seqs = new LinkedHashSet<>();
        for(int i=0; i<n; i++) {
            DNASequence mut_seq = mutate_seq(seq);
            mut_seq.setAccession(new AccessionID(String.format("%s_%d",seq_name,i)));
            out_seqs.add(mut_seq);
        }
        return out_seqs;
    }

    private DNASequence mutate_seq(DNASequence seq) throws CompoundNotFoundException {
        StringBuilder mut_seq = new StringBuilder();
        for (int i=0; i<seq.getLength(); i++) {
            Character nuc = seq.getCompoundAt(i+1).toString().charAt(0);
            if (Math.random() > nuc_mut_prob) {
                mut_seq.append(nuc);
                continue;
            }
            Character mut_nuc = nucs.get(rand.nextInt(nucs.size()-1));
            if (nuc == mut_nuc) nuc = nucs.get(3);
            mut_seq.append(nuc);
        }
        return new DNASequence(mut_seq.toString(), AmbiguityDNACompoundSet.getDNACompoundSet());
    }
}
