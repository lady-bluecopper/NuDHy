package samosa.test;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.jgrapht.graph.*;

import java.io.*;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import org.jgrapht.alg.isomorphism.VF2GraphIsomorphismInspector;
import samosa.structures.DirectedBipartiteGraph;
import samosa.utils.CMDLineParser;
import samosa.utils.Config;
import samosa.utils.io.Transformer;

/**
 * This class performs the graph isomorphism test for a list of random graphs 
 * generated using the sampler samplerType. The random graphs are read from
 * directoryPath. 
 * @author giulia
 */
public class IsomorphismTest {
    

    public static void main(String[] args) throws IOException {
        
        CMDLineParser.parse(args);
        
        // READ DATASET
        final String filePath = Path.of(Config.datasetsDir, Config.datasetName).toString();
        final Transformer T = new Transformer();
        final DirectedBipartiteGraph dirHyperGraph = samosa.utils.io.Reader.readDirectedHyperGraph(filePath, T);
        SimpleDirectedGraph<String, DefaultEdge> g = fromBipDirToSimplDir(dirHyperGraph);
        
        // LOAD SAMPLE FILE NAMES
        String dataName = Config.datasetName.split("\\.")[0];
        final String directoryPath = Path.of(Config.resultsDir, "samples", "NuDHy", dataName).toString();
        final File directory = new File(directoryPath);
        final File[] files = directory.listFiles();
        Map<String, List<String>> files_per_sampler = Maps.newHashMap();
        if (files != null) {
            for (File file : files) {
                if (file.isFile()) {
                    final String fileName = file.getName();
                    if (fileName.contains(dataName)) {
                        Map<String, String> fieldValueMap = samosa.utils.io.Reader.parseOutputFileName(fileName);
                        if (!fieldValueMap.get("sampleNumber").equals(String.valueOf(Config.seed))) {
                            continue;
                        }
                        String sampler = fieldValueMap.get("algorithm");
                        if (sampler.equals(Config.samplerType)) {
                            if (!files_per_sampler.containsKey(sampler)) {
                                files_per_sampler.put(sampler, Lists.newArrayList());
                            }
                            files_per_sampler.get(sampler).add(file.getAbsolutePath());
                        }
                    }
                }
            }
        }
        // RUN ISOMORPHISM
        for (Entry<String, List<String>> filesList : files_per_sampler.entrySet()) {
            System.out.println(filesList.getKey());
            runIsomorphism(g, filesList.getValue(), T);
        }
    }
    
    /**
     * Converts a DirBipGraph to a SimpleDirectedGraph.
     * @param dirHyperGraph directed bipartite graph
     * @return the graph as a SimpleDirectedGraph 
     */
    private static SimpleDirectedGraph fromBipDirToSimplDir(DirectedBipartiteGraph dirHyperGraph) {
        SimpleDirectedGraph<String, DefaultEdge> g = new SimpleDirectedGraph(DefaultEdge.class);
        // add the vertices
        for (int v = 0; v < dirHyperGraph.getLeftP().length; v++) {
            g.addVertex("L" + v);
        }
        for (int v = 0; v < dirHyperGraph.getRightP().length; v++) {
            g.addVertex("R" + v);
        }
        // add the edges
        for (int v = 0; v < dirHyperGraph.getLeftP().length; v++) {
            for (int u : dirHyperGraph.getLeftP()[v]) {
                g.addEdge("L" + v, "R" + u);
            }
        }
        for (int v = 0; v < dirHyperGraph.getLeftM().length; v++) {
            for (int u : dirHyperGraph.getLeftM()[v]) {
                g.addEdge("R" + u, "L" + v);
            }
        }
        return g;
    }
    
    /**
     * Run the isomorphism test between g and each graph in filePaths
     * @param g SimpleDirectedGraph to comapre with the graphs in filePaths
     * @param filePaths path to each sample graph
     * @param T transformer to get the outer node ids
     * @throws IOException 
     */
    private static void runIsomorphism(SimpleDirectedGraph g,
            List<String> filePaths,
            Transformer T) throws IOException {

        // initialize statistics
        Object2IntOpenHashMap<String> counts = new Object2IntOpenHashMap();
        counts.defaultReturnValue(0);
        List<String> samples = Lists.newArrayList();
        int counter = 0;
        DirectedBipartiteGraph dirHyperGraph;
        DirectedBipartiteGraph othG;
        SimpleDirectedGraph sample;
        SimpleDirectedGraph other;
        boolean visited;
        VF2GraphIsomorphismInspector inspector;
        for (String filePath : filePaths) {
            String[] lst = filePath.split("/");
            System.out.println(lst[lst.length-1]);
            counter += 1;
            // load hypergraph
            dirHyperGraph = samosa.utils.io.Reader.readDirectedHyperGraph(filePath, T);
            sample = fromBipDirToSimplDir(dirHyperGraph);
            visited = false;
            inspector = new VF2GraphIsomorphismInspector(sample, g);
            if (inspector.isomorphismExists()) {
                visited = true;
                counts.addTo("input", 1);
            } else {
                for (String othPath : samples) {
                    othG = samosa.utils.io.Reader.readDirectedHyperGraph(othPath, T);
                    other = fromBipDirToSimplDir(othG);
                    inspector = new VF2GraphIsomorphismInspector(sample, other);
                    if (inspector.isomorphismExists()) {
                        visited = true;
                        counts.addTo(othPath, 1);
                        break;
                    }
                }
            }
            if (!visited) {
                samples.add(filePath);
            }
            System.out.println(samples.size() + " different out of " + counter + " processed.");
        }
        System.out.println(Arrays.toString(counts.values().toIntArray()));
    }
}
