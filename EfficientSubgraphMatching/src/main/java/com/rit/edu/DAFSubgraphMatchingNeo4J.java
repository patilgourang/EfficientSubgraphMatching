package com.rit.edu;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.SimpleGraph;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

public class DAFSubgraphMatchingNeo4J {
	
	private static String queryFolderPath = "Proteins/query";
	private static String DBPath = "ProtiensDB";
	
	private static GraphDatabaseFactory dbFactory;
	private static GraphDatabaseService dbService;
	
	private static void initDB() {
		dbFactory = new GraphDatabaseFactory();
		dbService = dbFactory.newEmbeddedDatabase(new File(DBPath));
	}
	
	private static void shutdownDB() {
		dbService.shutdown();
	}
	public Graph<Vertex, DefaultEdge> createQueryGraph(String queryFileName) throws IOException {
		//SimpleGraph<Vertex, DefaultEdge> queryGraph = new SimpleGraph<Vertex, DefaultEdge>(DefaultEdge.class);
		Graph<Vertex, DefaultEdge> queryGraph = new SimpleGraph<>(DefaultEdge.class);
		File queryFile = new File(queryFolderPath+"/"+queryFileName);
		System.out.println("Query File Name : "+queryFileName);
		LinkedList<String> linesList = new LinkedList<>();
		BufferedReader bf = new BufferedReader(new FileReader(queryFile));
		String line = null;
		//int numOfNodes = Integer.parseInt(queryFileName.split("\\.")[1]);
		while ((line = bf.readLine()) != null)
        {
			linesList.add(line);
        }
		bf.close();
		
		int nodeCounter = Integer.parseInt(linesList.getFirst());
		linesList.removeFirst();
		Map<String, Vertex> vertexMap = new HashMap<String, Vertex>(); 
		for(int i = 0;i < nodeCounter; i++)
		{
			String labelArr [] = linesList.getFirst().split(" ");
			Vertex v = new Vertex("u" + labelArr[0], labelArr[1]);
			vertexMap.put("u" + labelArr[0], v);
			queryGraph.addVertex(v);
			linesList.removeFirst();
		}
		
		while(!linesList.isEmpty())
		{
			int relationCounter = Integer.parseInt(linesList.getFirst());
			linesList.removeFirst();
			for(int i = 0; i< relationCounter;i++)
			{
				String relationArr [] = linesList.getFirst().split(" ");
				queryGraph.addEdge(vertexMap.get("u" + relationArr[0]),vertexMap.get("u" + relationArr[1]));
				linesList.removeFirst();
			}
		}
		System.out.println("Query Graph"+ queryGraph);
		return queryGraph;
	}
	
	public void DAF(Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel)
	{
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = buildDAG(queryGraph, dataGraphLabel);
		buildCS(queryGraph, queryDAG, dataGraphLabel);
	}
	
	public SimpleDirectedGraph<Vertex, DefaultEdge> buildDAG(Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel)
	{
		Vertex root = selectQueryRoot(queryGraph, dataGraphLabel);
		//System.out.println(root);
		List<List<Vertex>> verticesByLevel = traverseGraphinBFSOrder(root, queryGraph);
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = buildDAGFromVertices(verticesByLevel, queryGraph, dataGraphLabel);
		System.out.println("Query DAG: " + queryDAG);
		return queryDAG;
	}
	
	public SimpleDirectedGraph<Vertex, DefaultEdge> reverseDAG(Graph<Vertex, DefaultEdge> queryGraph) {
		SimpleDirectedGraph<Vertex, DefaultEdge> reverseDAG = new SimpleDirectedGraph<>(DefaultEdge.class);
		for(Vertex v : queryGraph.vertexSet())
		{
			reverseDAG.addVertex(v);
		}
		for(Vertex v : queryGraph.vertexSet())
		{
			for(DefaultEdge edge : queryGraph.edgesOf(v)) {
				reverseDAG.addEdge(queryGraph.getEdgeTarget(edge), queryGraph.getEdgeSource(edge));
			}
		}
		
		System.out.println("Reverrse DAG" + reverseDAG);
		return reverseDAG;
	}
	
	private Vertex selectQueryRoot(Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel) {
		Vertex root = null;
		float value = Float.MAX_VALUE;
		for(Vertex queryVertex : queryGraph.vertexSet())
		{
			float currValue = (float) getInitialCandidateSetCount(queryVertex, queryGraph, dataGraphLabel) / queryGraph.edgesOf(queryVertex).size();
			if (currValue < value) {
				root = queryVertex;
				value = currValue;
			}
		}
		return root;
	}
	
	public int getInitialCandidateSetCount(Vertex queryVertex, Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel){
		int candidateSetCount = 0;
		String query = "MATCH(N1:" + dataGraphLabel + ":" + queryVertex.label +")-[:HASCONNECTION]-(N2) WITH COUNT(N2) AS NodeCount, N1 WHERE NodeCount >= "+ queryGraph.edgesOf(queryVertex).size() +" RETURN COUNT(N1);";
		Result rs = dbService.execute(query);
		while(rs.hasNext())
		{
			Map<String,Object> next = rs.next();
			for(Entry<String, Object> entry : next.entrySet())
			{
				//System.out.println(entry.getKey() +" : "+entry.getValue().toString());
				candidateSetCount = Integer.parseInt(entry.getValue().toString());
			}
		}
		rs.close();
		return candidateSetCount;
	}
	
	public Set<Integer> getInitialCandidateSet(Vertex queryVertex, Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel){
		Set<Integer> candidateSet = new LinkedHashSet<>();
		String query = "MATCH(N1:" + dataGraphLabel + ":" + queryVertex.label +")-[:HASCONNECTION]-(N2) WITH COUNT(N2) AS NodeCount, N1 WHERE NodeCount >= "+ queryGraph.edgesOf(queryVertex).size() +" RETURN N1.id;";
		Result rs = dbService.execute(query);
		while(rs.hasNext())
		{
			Map<String,Object> next = rs.next();
			for(Entry<String, Object> entry : next.entrySet())
			{
				candidateSet.add(Integer.parseInt(entry.getValue().toString()));
			}
		}
		rs.close();
		return candidateSet;
	}
	
	public Set<Integer> getAdjacentNodesOf(Integer node, Vertex queryVertex, String dataGraphLabel, String childLabel){
		Set<Integer> adjacentNodes = new TreeSet<>();
		String query = "MATCH(N1:" + dataGraphLabel + ":" + queryVertex.label +")-[:HASCONNECTION]-(N2:"+ childLabel +") WHERE N1.id = '"+ node +"' RETURN N2.id;";
		//System.out.println(query);
		Result rs = dbService.execute(query);
		while(rs.hasNext())
		{
			Map<String,Object> next = rs.next();
			for(Entry<String, Object> entry : next.entrySet())
			{
				adjacentNodes.add(Integer.parseInt(entry.getValue().toString()));
			}
		}
		rs.close();
		return adjacentNodes;
	}
	
	private List<List<Vertex>> traverseGraphinBFSOrder(Vertex startVertex, Graph<Vertex, DefaultEdge> queryGraph) {
		
		List<List<Vertex>> verticesByLevel = new ArrayList<>();
		Set<Vertex> visited = new HashSet<>(); 
		
		Queue<Vertex> queue = new LinkedList<>();
		queue.add(startVertex);
		visited.add(startVertex);
		
		while(!queue.isEmpty()) {
			List<Vertex> vertexList = new ArrayList<>();
			int size = queue.size();
			while(size != 0)
			{
				Vertex v = queue.poll();
				//System.out.print(v + " ");
				vertexList.add(v);
				for(DefaultEdge edge : queryGraph.edgesOf(v)) {
					Vertex target = queryGraph.getEdgeTarget(edge);
					if(target.equals(v)) {
						target = queryGraph.getEdgeSource(edge);
					}
					if (!visited.contains(target)) {
						queue.add(target);
						visited.add(target);
					}
				}
				size--;
			}
			verticesByLevel.add(vertexList);
			//System.out.println();
		}
		return verticesByLevel;
	}
	
	private SimpleDirectedGraph<Vertex, DefaultEdge> buildDAGFromVertices(List<List<Vertex>> verticesByLevel, Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel) {
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = new SimpleDirectedGraph<>(DefaultEdge.class);
		
		for(List<Vertex> vertexList : verticesByLevel) {
			for(Vertex v : vertexList) {
				queryDAG.addVertex(v);
			}
		}
		
		for(int i = 0; i < verticesByLevel.size() - 1; i++) {
			List<Vertex> upperLevel = verticesByLevel.get(i);
			List<Vertex> lowerLevel = verticesByLevel.get(i+1);
			
			//Direct edges from upper to lower levels
			for (int j = 0; j < upperLevel.size(); j++) {
				Vertex upperVertex = upperLevel.get(j);
				for (int k=0; k< lowerLevel.size(); k++) {
					Vertex lowerVertex = lowerLevel.get(k);
					if(queryGraph.containsEdge(upperVertex, lowerVertex) || queryGraph.containsEdge(lowerVertex, upperVertex)) {
						queryDAG.addEdge(upperVertex, lowerVertex);
					}
				}
			}
			
			// direct edges for vertices on same levels
			List<List<Vertex>> labelsGroup = groupAndSortVertices(lowerLevel, queryGraph, dataGraphLabel);
			int groups = labelsGroup.size();
			if(groups > 1) {
				for(int x=0; x < groups -1; x++) {
					List<Vertex> group1 = labelsGroup.get(x);
					List<Vertex> group2 = labelsGroup.get(x+1);
					for(int y = 0; y < group1.size(); y++) {
						Vertex v1 = group1.get(y);
						for(int z = 0; z < group2.size(); z++) {
							Vertex v2 = group2.get(z);
							if(queryGraph.containsEdge(v1, v2) || queryGraph.containsEdge(v2, v2)) {
								queryDAG.addEdge(v1, v2);
							}
						}
					}
				}
			}
			//System.out.println(labelsGroup);
		}
		//System.out.println("Query DAG" + queryDAG);
		return queryDAG;
	}
	
	/**
	 * Vertices by grouped by Labels. In each group the vertices are sorted in descending order of vertex degrees.  
	 * @param lowerLevel
	 * @param queryGraph
	 * @param dataGraph
	 */
	private List<List<Vertex>> groupAndSortVertices(List<Vertex> lowerLevel, Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel) {
		Map<String, List<Vertex>> labelsMap = new HashMap<>();
		for(Vertex v: lowerLevel)
		{
			if(labelsMap.containsKey(v.label))
			{
				List<Vertex> vertexList = labelsMap.get(v.label);
				vertexList.add(v);
				labelsMap.put(v.label, vertexList);
			} else {
				List<Vertex> vertexList = new ArrayList<>();
				vertexList.add(v);
				labelsMap.put(v.label, vertexList);
			}
		}
		
		List<List<Vertex>> labelsGroup = new ArrayList<>();
		List<String> labelOrder = getLabelOrder(dataGraphLabel, labelsMap.keySet());
		for(String label : labelOrder)
		{
			if(labelsMap.containsKey(label)) {
				List<Vertex> sortedList = sortVerticesByDegree(labelsMap.get(label), queryGraph);
				labelsGroup.add(sortedList);
			}
		}
		return labelsGroup;
	}
	
	/**
	 * Get the order of labels in the data graph (Most infrequent labels in graph come earlier)
	 * @param dataGraph
	 */
	private List<String> getLabelOrder(String dataGraphLabel, Set<String> labelsSet) {
		Map<String, Integer> labelCountsMap = new HashMap<>();
		int labelCount = 0;
		for(String label : labelsSet) {
			String query = "MATCH(N1:" + dataGraphLabel + ":" + label +") RETURN COUNT(N1);";
			Result rs = dbService.execute(query);
			while(rs.hasNext())
			{
				Map<String,Object> next = rs.next();
				for(Entry<String, Object> entry : next.entrySet())
				{
					//System.out.println(entry.getKey() +" : "+entry.getValue().toString());
					labelCount = Integer.parseInt(entry.getValue().toString());
				}
			}
			labelCountsMap.put(label, labelCount);
			rs.close();
		}
		
		LinkedHashMap<String, Integer> sorted = labelCountsMap
		        .entrySet()
		        .stream()
		        .sorted(Map.Entry.comparingByValue())
		        .collect(
		            Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2,
		                LinkedHashMap::new));
		
		return new ArrayList<>(sorted.keySet());
	}
	
	/**
	 * Sort vertices in descending order of their degree.
	 * @param verticesList
	 * @param queryGraph
	 * @return
	 */
	private List<Vertex> sortVerticesByDegree(List<Vertex> verticesList, Graph<Vertex, DefaultEdge> queryGraph) {
		verticesList.sort((Vertex v1, Vertex v2) -> queryGraph.edgesOf(v2).size() - queryGraph.edgesOf(v1).size());
		return verticesList;
	}
	
	private void buildCS(Graph<Vertex, DefaultEdge> queryGraph, SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, String dataGraphLabel)
	{
		SimpleDirectedGraph<Vertex, DefaultEdge> reverseQueryDAG = reverseDAG(queryGraph);
		
		Map<String, Set<Integer>> initialCS = new LinkedHashMap<>();
		for(Vertex vertex : queryDAG.vertexSet())
		{
			Set<Integer> candidateSet = getInitialCandidateSet(vertex, queryGraph, dataGraphLabel);
			initialCS.put(vertex.id, candidateSet);
			//System.out.println("Vertex: " + vertex + " Candidate Set: " + candidateSet);
		}
		System.out.println("Initial CS: " + initialCS);
		
		Map<String, Set<Integer>> refinedCS = DAGGraphDP(reverseQueryDAG, initialCS, dataGraphLabel);
		System.out.println("First Refinement" + refinedCS);
		refinedCS = DAGGraphDP(queryDAG, refinedCS, dataGraphLabel);
		System.out.println("Second Refinement" + refinedCS);
		refinedCS = DAGGraphDP(reverseQueryDAG, refinedCS, dataGraphLabel);
		System.out.println("Third Refinement" + refinedCS);
		refinedCS = DAGGraphDP(queryDAG, refinedCS, dataGraphLabel);
		System.out.println("Forth Refinement" + refinedCS);
	}
	
	private Map<String, Set<Integer>> DAGGraphDP(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Map<String, Set<Integer>> initialCS, String dataGraphLabel) {
		// TODO Auto-generated method stub
		Set<Vertex> visited = new HashSet<>();
		//Stack<Vertex> stack = new Stack<>();
		Queue<Vertex> queue = new LinkedList<>();
		Map<String, Set<Integer>> refinedCS = new HashMap<>();
		Vertex startVertex = findVertexWithInDegreeZero(queryDAG);
		//reverseTopologicalSort(queryDAG, visited, queue, queryDAG.vertexSet().stream().filter(v -> v.id.equals("u4")).findFirst().get());
		reverseTopologicalSort(queryDAG, visited, queue, startVertex);
		System.out.println("QueryDAG :" + queryDAG + "Start Vertex" + startVertex);
		//System.out.println(stack);
		/*while(!stack.empty()) {
			System.out.println(stack.pop());
		}*/
		while (!queue.isEmpty()) {
			//System.out.println(queue.poll());
			Vertex v = queue.poll();
			Set<DefaultEdge> outgoingEdges = queryDAG.outgoingEdgesOf(v);
			//System.out.println("Vertex: " + v + " Outgoing Edges :" + outgoingEdges);
			Set<Vertex> children = outgoingEdges.stream().map(edge -> queryDAG.getEdgeTarget(edge)).collect(Collectors.toSet());
			System.out.println("Vertex: " + v + " Children :" + children);
			if (children.isEmpty()) {
				refinedCS.put(v.id, initialCS.get(v.id));
			}
			else {
				for(Integer node : initialCS.get(v.id)) {
					boolean flag = true;
					for(Vertex child : children) {
						Set<Integer> neighbours = getAdjacentNodesOf(node, v, dataGraphLabel, child.label);
						
						
						System.out.println("Node:" + node + " -> Child:" + child + ":" + neighbours);
						/*if(!refinedCS.get(child.id).containsAll(neighbours)) {
							flag = false;
						}*/
						neighbours.retainAll(refinedCS.get(child.id));
						if(neighbours.isEmpty()) {
							flag = false;
						}
						
					}
					if(flag) {
						if(refinedCS.containsKey(v.id)) {
							Set<Integer> candidates = refinedCS.get(v.id);
							candidates.add(node);
							refinedCS.put(v.id, candidates);
						}
						else {
							Set<Integer> candidates = new TreeSet<>();
							candidates.add(node);
							refinedCS.put(v.id, candidates);
						}
					}
					else {
						System.out.println("Refined CS: " + refinedCS);
						System.out.println("Removing " + node);
					}
				}
			}
		}
		System.out.println("Refined CS: " + refinedCS);
		return refinedCS;
	}

	private Vertex findVertexWithInDegreeZero(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG) {
		for(Vertex v: queryDAG.vertexSet())
		{
			if(queryDAG.incomingEdgesOf(v).size() == 0)
				return v;
		}
		throw new IllegalStateException("no vertex with In Degree 0 found");
	}

	private void reverseTopologicalSort(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Set<Vertex> visited, Queue<Vertex> queue, Vertex v) {
		visited.add(v);
		
		for(DefaultEdge edge:queryDAG.edgesOf(v))
		{
			Vertex target = queryDAG.getEdgeTarget(edge);
			if(target.equals(v)) {
				//target = queryDAG.getEdgeSource(edge);
				continue;
			}
			if (!visited.contains(target)) {
				reverseTopologicalSort(queryDAG, visited, queue, target);
			}
		}
		queue.offer(v);
		
	}
	
	private void topologicalSort(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Set<Vertex> visited, Stack<Vertex> stack, Vertex v) {
		// TODO Auto-generated method stub
		visited.add(v);
		
		for(DefaultEdge edge:queryDAG.edgesOf(v))
		{
			Vertex target = queryDAG.getEdgeTarget(edge);
			if(target.equals(v)) {
				target = queryDAG.getEdgeSource(edge);
			}
			if (!visited.contains(target)) {
				topologicalSort(queryDAG, visited, stack, target);
			}
		}
		stack.push(v);
		
	}

	public static void main(String[] args) throws IOException {
		initDB();
		DAFSubgraphMatchingNeo4J daf = new DAFSubgraphMatchingNeo4J();
		Graph<Vertex, DefaultEdge> queryGraph = daf.createQueryGraph("a_testquery.4.sub.grf");
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = daf.buildDAG(queryGraph, "a_test");
		daf.DAF(queryGraph, "a_test");
		/*SimpleDirectedGraph<Vertex, DefaultEdge> reverseQueryDAG = daf.reverseDAG(queryGraph);
		daf.buildCS(queryGraph, reverseQueryDAG, "a_test");*/
		shutdownDB();
	}
}
