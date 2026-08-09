package core;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import preprocessing.GeneTrees;
import taxon.Taxon;
import tree.RangeBipartition;
import tree.Tree;
import utils.Config;

class InferenceRegressionTest {

    @BeforeEach
    void resetGlobalTaxonIds() {
        Taxon.count = 0;
    }

    @Test
    void smallFixturePreservesBaselineScoreAndTopology() throws Exception {
        Config.ComputationMode oldMode = Config.COMPUTATION_MODE;
        try {
            Config.COMPUTATION_MODE = Config.ComputationMode.CPU_SINGLE;
            Path input = Path.of("tests", "resources", "regression-small.tre");
            GeneTrees geneTrees = new GeneTrees(input.toString());
            geneTrees.readTaxaNames();
            geneTrees.readGeneTrees(null);
            List<RangeBipartition> candidates = geneTrees.generateExtendedCandidateBipartitions(false);

            InferenceDP inference = new InferenceDP(geneTrees, candidates);
            assertEquals(47.0, inference.solve());

            Tree expected = new Tree("((D,E),(C,(A,B)));", geneTrees.taxaMap);
            assertEquals(true, expected.equalsStructure(inference.reconstructTree()));
            assertEquals(13, geneTrees.rangeBipartitions.size());
        } finally {
            Config.COMPUTATION_MODE = oldMode;
        }
    }

    @Test
    void adaptiveBitsetPathPreservesBaselineScoreAndTopology() throws Exception {
        Config.ComputationMode oldMode = Config.COMPUTATION_MODE;
        try {
            Config.COMPUTATION_MODE = Config.ComputationMode.CPU_SINGLE;
            Path input = Path.of("tests", "resources", "regression-bitset.tre");
            GeneTrees geneTrees = new GeneTrees(input.toString());
            geneTrees.readTaxaNames();
            geneTrees.readGeneTrees(null);
            List<RangeBipartition> candidates = geneTrees.generateExtendedCandidateBipartitions(false);

            InferenceDP inference = new InferenceDP(geneTrees, candidates);
            assertEquals(433.0, inference.solve());

            Tree expected = new Tree("(H,(G,((E,F),((C,D),(A,B)))));", geneTrees.taxaMap);
            assertEquals(true, expected.equalsStructure(inference.reconstructTree()));
            assertEquals(29, geneTrees.rangeBipartitions.size());
        } finally {
            Config.COMPUTATION_MODE = oldMode;
        }
    }
}
