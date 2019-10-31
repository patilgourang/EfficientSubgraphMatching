package com.rit.edu;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.jgrapht.graph.SimpleGraph;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import static java.util.stream.Collectors.*;
import static java.util.Map.Entry.*;

public class DAFSubgraphMatchingNeo4J {
	
	private static String queryFolderPath = "Proteins/query";
	private static String DBPath = "ProtiensDB";
	private static String groundTruthFolderPath = "Proteins/ground_truth";
	private static String dataFolderPath = "Proteins/target";
	
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
		//System.out.println("Query File Name : "+queryFileName);
		
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
		//System.out.println("Query Graph: "+ queryGraph);
		return queryGraph;
	}
	
	public List<String> DAF(Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel)
	{
		HashMap<Integer,ArrayList<Integer>> neighboursMap = createNeighboursMap(dataGraphLabel);
		//System.out.println("Data Graph Name : " + dataGraphLabel);
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = buildDAG(queryGraph, dataGraphLabel);
		List<Vertex> DAGorder = topologicalSort(queryDAG);
		Map<String, Set<Integer>> CS = buildCS(queryGraph, queryDAG, dataGraphLabel, neighboursMap);
		long startTimeBuildCS = System.currentTimeMillis();//System.nanoTime();
		Map<Integer, Map<String, Set<Integer>>> adjacencyList = buildCSAdjacencyList(CS, queryDAG, dataGraphLabel, neighboursMap);
		long buildCSTime = (System.currentTimeMillis()-startTimeBuildCS);
		System.out.println("Time taken for Bulding CS in ms :"+buildCSTime);
		long startTime = System.currentTimeMillis();//System.nanoTime();
		List<String> results = backtrack(new LinkedHashMap<>(), queryDAG, CS, DAGorder, new LinkedHashSet<>(), adjacencyList);
		long queryTime = (System.currentTimeMillis()-startTime);
		System.out.println("Time taken for Backtracking in ms :"+queryTime);
		//System.out.println(results);
		return results;
	}
	
	public List<String> DAFFailingSets(Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel)
	{
		HashMap<Integer,ArrayList<Integer>> neighboursMap = createNeighboursMap(dataGraphLabel);
		//System.out.println("Data Graph Name : " + dataGraphLabel);
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = buildDAG(queryGraph, dataGraphLabel);
		List<Vertex> DAGorder = topologicalSort(queryDAG);
		Map<String, Set<Integer>> CS = buildCS(queryGraph, queryDAG, dataGraphLabel, neighboursMap);
		long startTimeBuildCS = System.currentTimeMillis();//System.nanoTime();
		Map<Integer, Map<String, Set<Integer>>> adjacencyList = buildCSAdjacencyList(CS, queryDAG, dataGraphLabel, neighboursMap);
		long buildCSTime = (System.currentTimeMillis()-startTimeBuildCS);
		System.out.println("Time taken for Bulding CS in ms :"+buildCSTime);
		long startTime = System.currentTimeMillis();//System.nanoTime();
		Map<String, Set<String>> ancestorsmap = getAncestors(queryDAG);
		List<String> results = backtrackFailingSets(new LinkedHashMap<>(), queryDAG, CS, DAGorder, new LinkedHashSet<>(), adjacencyList, ancestorsmap, new HashMap<>(), new HashMap<>());
		long queryTime = (System.currentTimeMillis()-startTime);
		System.out.println("Time taken for Backtracking in ms :"+queryTime);
		//System.out.println(results);
		return results;
	}
	
	private Map<Integer, Map<String, Set<Integer>>> buildCSAdjacencyList(Map<String, Set<Integer>> CS, SimpleDirectedGraph<Vertex,DefaultEdge> queryDAG, String dataGraphLabel, HashMap<Integer,ArrayList<Integer>> neighboursMap) {
		// TODO Auto-generated method stub
		Map<Integer, Map<String, Set<Integer>>> adjacencyList = new LinkedHashMap<>();
		//Map<String, Set<Integer>> adjacencyList = new LinkedHashMap<>();
		for(Entry<String, Set<Integer>> entry : CS.entrySet())
		{
			Set<Integer> candidates_u = entry.getValue();
			for(int v : candidates_u) {
				Map<String, Set<Integer>> adjacencyMap = new LinkedHashMap<>();
				Vertex vertex = queryDAG.vertexSet().stream().filter(V -> V.id.equals(entry.getKey())).findFirst().orElse(null);
				//Set<Integer> adjacentNodesOfV = getAdjacentNodes(v, vertex, dataGraphLabel); /// here
				Set<Integer> adjacentNodesOfV = new HashSet<Integer>(neighboursMap.get(v));
				for(DefaultEdge edge : queryDAG.outgoingEdgesOf(vertex)) {
					Vertex uc = queryDAG.getEdgeTarget(edge);
					Set<Integer> result = adjacentNodesOfV.stream()
							  .distinct()
							  .filter(CS.get(uc.id)::contains)
							  .collect(Collectors.toSet());
					//System.out.println();
					adjacencyMap.put(vertex.id + "_" + uc.id, result);
				}
				if(!adjacencyMap.isEmpty()) {
					if(adjacencyList.containsKey(v)) {
						Map<String, Set<Integer>> map = adjacencyList.get(v);
						map.putAll(adjacencyMap);
						adjacencyList.put(v, map);
					} else {
						adjacencyList.put(v, adjacencyMap);
					}
				}
			}
		}
		//System.out.println(adjacencyList);
		return adjacencyList;
	}

	public SimpleDirectedGraph<Vertex, DefaultEdge> buildDAG(Graph<Vertex, DefaultEdge> queryGraph, String dataGraphLabel)
	{
		Vertex root = selectQueryRoot(queryGraph, dataGraphLabel);
		//System.out.println(root);
		List<List<Vertex>> verticesByLevel = traverseGraphinBFSOrder(root, queryGraph);
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = buildDAGFromVertices(verticesByLevel, queryGraph, dataGraphLabel);
		//System.out.println("Query DAG: " + queryDAG);
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
		
		//System.out.println("Reverse DAG: " + reverseDAG);
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
		//System.out.println("Initial Candidate Set" + queryVertex.label + " : " + candidateSet );
		return candidateSet;
	}
	
	public static HashMap<Integer,ArrayList<Integer>> createNeighboursMap(String fileLabel)
	{
		LinkedHashMap<Integer,ArrayList<Integer>> neighboursMap = new LinkedHashMap<>();
		//String query = "MATCH(N1:"+fileLabel+") RETURN N1,N1.neighbours";
		String query = "MATCH(N1:" + fileLabel +")-[:HASCONNECTION]-(N2) RETURN N1, collect(N2.id);";
		//System.out.println(query);
		Result rs = dbService.execute(query);
		
		while(rs.hasNext())
		{
			int nodeId = -1;
			ArrayList<Integer> neighbourList = new ArrayList<>();
			Map<String,Object> next = rs.next();
			for(Entry<String, Object> entry : next.entrySet())
			{
				if(entry.getKey().equals("N1"))
				{
					String node = entry.getValue().toString();
					nodeId = Integer.parseInt(node.substring(5, node.length()-1)) % 10000;
				}
				else
				{
					String neighbours = entry.getValue().toString();
					neighbours = neighbours.replaceAll("\\[", "").replaceAll("\\]","");
					neighbours = neighbours.replaceAll("\\s", "");
					String arr [] = neighbours.split(",");
					for(int i=0;i<arr.length;i++)	
					{
						neighbourList.add(Integer.parseInt(arr[i]));
					}
				}
			}
			neighboursMap.put(nodeId, neighbourList);
		}
		rs.close();
		return neighboursMap;
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
	
	public Set<Integer> getAdjacentNodes(Integer node, Vertex queryVertex, String dataGraphLabel){
		Set<Integer> adjacentNodes = new TreeSet<>();
		String query = "MATCH(N1:" + dataGraphLabel +")-[:HASCONNECTION]-(N2) WHERE N1.id = '"+ node +"' RETURN N2.id;";
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
	
	private Map<String, Set<Integer>> buildCS(Graph<Vertex, DefaultEdge> queryGraph, SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, String dataGraphLabel, HashMap<Integer,ArrayList<Integer>> neighboursMap)
	{
		
		SimpleDirectedGraph<Vertex, DefaultEdge> reverseQueryDAG = reverseDAG(queryDAG);
		
		Map<String, Set<Integer>> initialCS = new TreeMap<>();
		for(Vertex vertex : queryDAG.vertexSet())
		{
			Set<Integer> candidateSet = getInitialCandidateSet(vertex, queryGraph, dataGraphLabel);
			initialCS.put(vertex.id, candidateSet);
			//System.out.println("Vertex: " + vertex + " Candidate Set: " + candidateSet);
		}
		//System.out.println();
		//System.out.println("Initial CS: " + initialCS);
		
		//System.out.print("Initial CS: ");
		StringBuilder initialCSStr = new StringBuilder("Initial CS: ");
		for(Entry<String, Set<Integer>> e : initialCS.entrySet()){
			//System.out.print(e.getKey() + " -> Candidate Set Size : (" + e.getValue().size() + ") ");
			//System.out.print(e.getKey() + ":" + e.getValue().size() + ", ");
			initialCSStr.append(e.getKey() + ":" + e.getValue().size() + ", ");
		}
		initialCSStr.delete(initialCSStr.length()-2,initialCSStr.length());
		System.out.println(initialCSStr);
		long startTime = System.currentTimeMillis();
		//System.out.println();
		Map<String, Set<Integer>> refinedCS = DAGGraphDP(reverseQueryDAG, initialCS, dataGraphLabel, neighboursMap);
		//System.out.println("First Refinement" + refinedCS);
		/*for(Entry<String, Set<Integer>> e : refinedCS.entrySet()){
			System.out.print(e.getKey() + ": Candidate Set Size : (" + e.getValue().size() + ") ; ");
		}*/
		//System.out.println();
		refinedCS = DAGGraphDP(queryDAG, refinedCS, dataGraphLabel, neighboursMap);
		//System.out.println("Second Refinement" + refinedCS);
		/*for(Entry<String, Set<Integer>> e : refinedCS.entrySet()){
			System.out.print(e.getKey() + ": Candidate Set Size : (" + e.getValue().size() + ") ; ");
		}*/
		//System.out.println();
		refinedCS = DAGGraphDP(reverseQueryDAG, refinedCS, dataGraphLabel, neighboursMap);
		//System.out.println("Third Refinement" + refinedCS);
		/*for(Entry<String, Set<Integer>> e : refinedCS.entrySet()){
			System.out.print(e.getKey() + ": Candidate Set Size : (" + e.getValue().size() + ") ; ");
		}*/
		//System.out.println();
		//refinedCS = DAGGraphDP(queryDAG, refinedCS, dataGraphLabel);
		//System.out.println("Forth Refinement" + refinedCS);
		/*for(Entry<String, Set<Integer>> e : refinedCS.entrySet()){
			System.out.print(e.getKey() + ": Candidate Set Size=(" + e.getValue().size() + ") ; ");
		}*/
		//System.out.println();
		
		long queryTime = (System.currentTimeMillis()-startTime);
		
		//System.out.print("Refined CS: ");
		StringBuilder refinedCSStr = new StringBuilder("Refined CS: ");
		for(Entry<String, Set<Integer>> e : refinedCS.entrySet()){
			//System.out.println(e.getKey() + " : " + e.getValue());
			//System.out.println(e.getKey() + " -> Size :" + e.getValue().size());
			//System.out.println(e.getKey() + " -> Candidate Set Size : (" + e.getValue().size() + ") ");
			//System.out.print(e.getKey() + ":" + e.getValue().size() + ", ");
			refinedCSStr.append(e.getKey() + ":" + e.getValue().size() + ", ");
		}
		refinedCSStr.delete(refinedCSStr.length()-2,refinedCSStr.length());
		System.out.println(refinedCSStr);
		System.out.println("Time taken for CS refinement in ms :"+queryTime);
		return refinedCS;
	}
	
	private Map<String, Set<Integer>> DAGGraphDP(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Map<String, Set<Integer>> initialCS, String dataGraphLabel, HashMap<Integer,ArrayList<Integer>> neighboursMap) {
		
		//HashMap<Integer,ArrayList<Integer>> neighboursMap = createNeighboursMap(dataGraphLabel);
		Map<String, Set<Integer>> refinedCS = new TreeMap();
		//Vertex startVertex = findVertexWithInDegreeZero(queryDAG);
		//reverseTopologicalSort(queryDAG, visited, queue, queryDAG.vertexSet().stream().filter(v -> v.id.equals("u4")).findFirst().get());
		//reverseTopologicalSortHelper(queryDAG, visited, queue, startVertex);
		//System.out.println("QueryDAG :" + queryDAG + "Start Vertex" + startVertex);
		//System.out.println(stack);
		/*while(!stack.empty()) {
			System.out.println(stack.pop());
		}*/
		Queue<Vertex> queue = reverseTopologicalSort(queryDAG);
		//System.out.println("Reverse Topological Order: "+queue);
		while (!queue.isEmpty()) {
			//System.out.println(queue.poll());
			Vertex v = queue.poll();
			Set<DefaultEdge> outgoingEdges = queryDAG.outgoingEdgesOf(v);
			//System.out.println("Vertex: " + v + " Outgoing Edges :" + outgoingEdges);
			Set<Vertex> children = outgoingEdges.stream().map(edge -> queryDAG.getEdgeTarget(edge)).collect(Collectors.toSet());
			//System.out.println("Vertex: " + v + " Children :" + children);
			if (children.isEmpty()) {
				refinedCS.put(v.id, initialCS.get(v.id));
			}
			else {
				for(Integer node : initialCS.get(v.id)) {
					boolean flag = true;
					for(Vertex child : children) {
						//Set<Integer> neighbours = getAdjacentNodesOf(node, v, dataGraphLabel, child.label);
						//System.out.println("Node:" + node + " -> Child:" + child + ":" + neighbours);
						/*if(!refinedCS.get(child.id).containsAll(neighbours)) {
							flag = false;
						}*/
						/*neighbours.retainAll(refinedCS.get(child.id));
						if(neighbours.isEmpty()) {
							flag = false;
						}*/
						//Getting neighbors from DB
						//Set<Integer> neighbours = getAdjacentNodesOf(node, v, dataGraphLabel, child.label);
						//System.out.println("Node:" + node + " -> Neighbours :" + neighbours);
						Set<Integer> neighbours = new HashSet<Integer>(neighboursMap.get(node)); 
						Set<Integer> candidateSetOfChild = refinedCS.get(child.id);
						//System.out.println("Child:" + child +" -> Candidate Set:" +candidateSetOfChild);
						neighbours.retainAll(candidateSetOfChild);
						if(neighbours.isEmpty()) {
							flag = false;
						}
						
						/*
						 * Approach Edge checking
						boolean childFlag = false;
						Set<Integer> candidateSetOfChild = refinedCS.get(child.id);
						if(candidateSetOfChild == null) {
							System.out.println("Child:" + child +" -> Candidate Set:" +candidateSetOfChild);
							System.out.println(refinedCS.keySet());
						}
							
						for(Integer candidate : candidateSetOfChild)
						{
							if(checkIfEdgeExist(node, candidate, dataGraphLabel))
							{
								childFlag = true;
								break;
							}
						}
						if(!childFlag) {
							flag = false;
							break;
						}*/
						
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
						//System.out.println("Adding Node" + node + "to" + v.id);
					}
					else {
						//System.out.println("Refined CS: " + refinedCS);
						//System.out.println("Removing " + node);
					}
				}
			}
			/*System.out.println("Refined CS: " + refinedCS.keySet());
			for(Entry<String, Set<Integer>> e : refinedCS.entrySet()){
				System.out.print(e.getKey() + ": Candidate Set Size : (" + e.getValue().size() + ") ; ");
			}
			System.out.println();*/
		}
		//System.out.println("Refined CS: " + refinedCS);
		return refinedCS;
	}

	private boolean checkIfEdgeExist(Integer node, Integer candidate, String dataGraphLabel) {
		String query = "RETURN EXISTS( (:"+ dataGraphLabel + " {id: '" + node + "'})-[:HASCONNECTION]-(:" + dataGraphLabel +" {id: '" + candidate + "'}) )";
		boolean flag = false;
		Result rs = dbService.execute(query);
		while(rs.hasNext())
		{
			Map<String,Object> next = rs.next();
			for(Entry<String, Object> entry : next.entrySet())
			{
				if(entry.getValue().toString().equals("true"))
					flag = true;
			}
		}
		rs.close();
		return flag;
	}

	private Vertex findVertexWithInDegreeZero(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG) {
		for(Vertex v: queryDAG.vertexSet())
		{
			if(queryDAG.incomingEdgesOf(v).size() == 0)
				return v;
		}
		throw new IllegalStateException("no vertex with In Degree 0 found");
	}
	
	private Queue<Vertex> reverseTopologicalSort(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG) 
    { 
		Set<Vertex> visited = new HashSet<>();
		Queue<Vertex> queue = new LinkedList<>(); 
  
        for (Vertex v : queryDAG.vertexSet()) 
            if (!visited.contains(v)) 
            	reverseTopologicalSortHelper(v, queryDAG, visited,  queue);
        return queue;
    }
	
	private void reverseTopologicalSortHelper(Vertex v, SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Set<Vertex> visited, Queue<Vertex> queue) {
		visited.add(v);
		
		for(DefaultEdge edge:queryDAG.edgesOf(v))
		{
			Vertex target = queryDAG.getEdgeTarget(edge);
			if(target.equals(v)) {
				//target = queryDAG.getEdgeSource(edge);
				continue;
			}
			if (!visited.contains(target)) {
				reverseTopologicalSortHelper(target, queryDAG, visited, queue);
			}
		}
		queue.offer(v);
	}
	
	private List<Vertex> topologicalSort(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG) 
    { 
		Set<Vertex> visited = new HashSet<>();
		Stack<Vertex> stack = new Stack<>(); 
  
        for (Vertex v : queryDAG.vertexSet()) 
            if (!visited.contains(v)) 
            	topologicalSortHelper(queryDAG, visited,  stack, v);
        
        List<Vertex> topologicalOrder = new ArrayList<>();
        while(!stack.isEmpty()) {
        	topologicalOrder.add(stack.pop());
        }
        return topologicalOrder;
    }
	
	private void topologicalSortHelper(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Set<Vertex> visited, Stack<Vertex> stack, Vertex v) {
		// TODO Auto-generated method stub
		visited.add(v);
		
		for(DefaultEdge edge:queryDAG.edgesOf(v))
		{
			Vertex target = queryDAG.getEdgeTarget(edge);
			if(target.equals(v)) {
				target = queryDAG.getEdgeSource(edge);
			}
			if (!visited.contains(target)) {
				topologicalSortHelper(queryDAG, visited, stack, target);
			}
		}
		stack.push(v);
		
	}
	
	/*private void getDAGOrdering(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG) {
		List<Vertex> order = topologicalSort(queryDAG);
		
	}*/
	
	private Set<Integer> getExtendableCandidates(String u_id, SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Map<Integer, Map<String, Set<Integer>>> adjacencyList, Map<String, Integer> partialEmbedding){
		Set<Integer> extendableCandidates = new LinkedHashSet<>();
		for(DefaultEdge edge : queryDAG.incomingEdgesOf(queryDAG.vertexSet().stream().filter(u -> u.id.equals(u_id)).findFirst().get())) {
			String parent = queryDAG.getEdgeSource(edge).id;
			Set<Integer> candidates =  adjacencyList.get(partialEmbedding.get(parent)).get(parent + "_" + u_id);
			if(extendableCandidates.isEmpty()) {
				if(candidates != null && !candidates.isEmpty()) {
					//System.out.println();
					extendableCandidates.addAll(candidates);
				}
				
			}
			else {
				extendableCandidates.retainAll(candidates);
			}
		}
		return extendableCandidates;
	}
	
	private List<String> backtrack(Map<String, Integer> partialEmbedding, SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Map<String, Set<Integer>> CS, List<Vertex> DAGorder, Set<Integer> visited, Map<Integer, Map<String, Set<Integer>>> adjacencyList)  {
		//System.out.println("Embedding" +partialEmbedding);
		List<String> results = new ArrayList<>();
		if(partialEmbedding.size() == queryDAG.vertexSet().size()) {
			Map<String, Integer> sorted = partialEmbedding
					.entrySet()
					.stream()
					.sorted(comparingByKey())
					.collect(
					toMap(Map.Entry::getKey, Map.Entry::getValue,
					(e1, e2) -> e2, LinkedHashMap::new));
			//System.out.println("Embedding" +sorted);
			results.add(convertEmbeddingToString(sorted));
			//System.out.println(convertEmbeddingToString(sorted));
		}
		else if(partialEmbedding.isEmpty() || partialEmbedding.size() == 0)
		{
			String nodeToProcess  = DAGorder.get(partialEmbedding.size()).id;
			for(int v : CS.get(nodeToProcess)){
				partialEmbedding.put(nodeToProcess, v);
				visited.add(v);
				results.addAll(backtrack(partialEmbedding, queryDAG, CS, DAGorder, visited, adjacencyList));
				partialEmbedding.remove(nodeToProcess);
				visited.remove((Integer)v);
			}
		}
		else 
		{
			String nodeToProcess  = DAGorder.get(partialEmbedding.size()).id;
			Set<Integer> extendableCandidates = getExtendableCandidates(nodeToProcess, queryDAG, adjacencyList, partialEmbedding);
			/*if(extendableCandidates.isEmpty()) {
				System.out.println("Embedding" + partialEmbedding + " emptyset-class->" + nodeToProcess);
			}*/
			for(int v : extendableCandidates)
			//for(int v : CS.get(nodeToProcess))
			{
				if(!visited.contains(v)) {
					partialEmbedding.put(nodeToProcess, v);
					visited.add(v);
					results.addAll(backtrack(partialEmbedding, queryDAG, CS, DAGorder, visited, adjacencyList));
					partialEmbedding.remove(nodeToProcess);
					visited.remove(v);
				} /*else {
					System.out.println("Embedding" + partialEmbedding + " conflict-class->" + nodeToProcess + ":" + v);
				}*/
			}
		}
		return results;
	}
	
	public HashMap<String, Set<String>> getAncestors(SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG) {
		HashMap<String, Set<String>> ancestorsMap = new HashMap<>();
		for(Vertex vertex : queryDAG.vertexSet()) {
			Set<String> ancestors = new HashSet<>();
			for(DefaultEdge edge : queryDAG.incomingEdgesOf(vertex)) {
				Vertex v = queryDAG.getEdgeSource(edge);
				ancestors.add(v.id);
			}
			ancestorsMap.put(vertex.id, ancestors);
		}
		return ancestorsMap;
	} 
	
	private List<String> backtrackFailingSets(Map<String, Integer> partialEmbedding, SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG, Map<String, Set<Integer>> CS, List<Vertex> DAGorder, Set<Integer> visited, Map<Integer, Map<String, Set<Integer>>> adjacencyList, Map<String, Set<String>> ancestorMap, Map<String, Set<String>> failingSets, Map<String, Set<String>> nodeChildrenMap)  {
		//System.out.println("Embedding" +partialEmbedding);
		List<String> results = new ArrayList<>();
		if(partialEmbedding.size() == queryDAG.vertexSet().size()) {
			Map<String, Integer> sorted = partialEmbedding
					.entrySet()
					.stream()
					.sorted(comparingByKey())
					.collect(
					toMap(Map.Entry::getKey, Map.Entry::getValue,
					(e1, e2) -> e2, LinkedHashMap::new));
			//System.out.println("Embedding" +sorted);
			results.add(convertEmbeddingToString(sorted));
			//System.out.println(convertEmbeddingToString(sorted));
		}
		else if(partialEmbedding.isEmpty() || partialEmbedding.size() == 0)
		{
			String nodeToProcess  = DAGorder.get(partialEmbedding.size()).id;
			for(int v : CS.get(nodeToProcess)){
				partialEmbedding.put(nodeToProcess, v);
				visited.add(v);
				results.addAll(backtrackFailingSets(partialEmbedding, queryDAG, CS, DAGorder, visited, adjacencyList, ancestorMap, failingSets, nodeChildrenMap));
				partialEmbedding.remove(nodeToProcess);
				visited.remove((Integer)v);
			}
		}
		else 
		{
			String nodeToProcess  = DAGorder.get(partialEmbedding.size()).id;
			Set<Integer> extendableCandidates = getExtendableCandidates(nodeToProcess, queryDAG, adjacencyList, partialEmbedding);
			if(extendableCandidates.isEmpty()) {
				//System.out.println("Embedding" + partialEmbedding + " emptyset-class->" + nodeToProcess);
			}
			int i = 1;
			for(int v : extendableCandidates)
			//for(int v : CS.get(nodeToProcess))
			{
				if(i > 1) {
					if(failingSets.containsKey(partialEmbedding.toString()))
						System.out.println(failingSets.get(partialEmbedding.toString()) + ":" + nodeToProcess);
				}
				if(!visited.contains(v)) {
					partialEmbedding.put(nodeToProcess, v);
					visited.add(v);
					results.addAll(backtrackFailingSets(partialEmbedding, queryDAG, CS, DAGorder, visited, adjacencyList, ancestorMap, failingSets, nodeChildrenMap));
					if(nodeChildrenMap.containsKey(partialEmbedding.toString())) {
						String childNode = checkChildNode(failingSets, nodeChildrenMap, partialEmbedding, DAGorder.get(partialEmbedding.size()).id);
						Set<String> failingSet = new HashSet<>();
						if(!childNode.equals("")) {
							failingSet = failingSets.get(childNode);
						} else {
							for(String child : nodeChildrenMap.get(partialEmbedding.toString())) {
								failingSet.addAll(failingSets.get(child));
							}
						}
						
						failingSets.put(partialEmbedding.toString(), failingSet);
						
						Map<String, Integer> parentPartialEmbedding = new LinkedHashMap<>(partialEmbedding);
						parentPartialEmbedding.remove(nodeToProcess);
						
						
						if(nodeChildrenMap.containsKey(parentPartialEmbedding.toString())) {
							Set<String> childrens = nodeChildrenMap.get(parentPartialEmbedding.toString());
							childrens.add(partialEmbedding.toString());
							nodeChildrenMap.put(parentPartialEmbedding.toString(), childrens);
						}
						else {
							Set<String> childrens = new HashSet<>();
							childrens.add(partialEmbedding.toString());
							nodeChildrenMap.put(parentPartialEmbedding.toString(), childrens);
						}
						
						/*if(failingSets.containsKey(parentPartialEmbedding.toString())) {
							Map<String, Set<String>> failingSetMap = failingSets.get(parentPartialEmbedding.toString());
							failingSetMap.put(parentPartialEmbedding.toString(), failingSet);
							failingSets.put(partialEmbedding.toString(), failingSetMap);
						} else {
							Map<String, Set<String>> failingSetMap = new LinkedHashMap<>();
							failingSetMap.put(parentPartialEmbedding.toString(), failingSet);
							failingSets.put(partialEmbedding.toString(), failingSetMap);
						}*/
						
						if(!failingSets.get(partialEmbedding.toString()).contains(nodeToProcess)) {
							//System.out.println("prune");
							partialEmbedding.remove(nodeToProcess);
							visited.remove(v);
							return results;
						}
						
					}
					partialEmbedding.remove(nodeToProcess);
					visited.remove(v);
				} else {
					//System.out.println("Embedding" + partialEmbedding + " conflict-class->" + nodeToProcess + ":" + v);
					String _u = partialEmbedding.entrySet().stream().filter(e -> e.getValue() == v).map(Map.Entry::getKey).findFirst().orElse(null);
					//System.out.println(_u + ":" + nodeToProcess);
					Map<String, Integer> failingSetNode = new LinkedHashMap<>(partialEmbedding);
					failingSetNode.put(nodeToProcess, v);
					Set<String> failingSet = new HashSet<>();
					failingSet.add(_u);
					failingSet.add(nodeToProcess);
					failingSet.addAll(ancestorMap.get(_u));
					failingSet.addAll(ancestorMap.get(nodeToProcess));
					
					failingSets.put(failingSetNode.toString(), failingSet);
					
					if(nodeChildrenMap.containsKey(partialEmbedding.toString())) {
						Set<String> childrens = nodeChildrenMap.get(partialEmbedding.toString());
						childrens.add(failingSetNode.toString());
						nodeChildrenMap.put(partialEmbedding.toString(), childrens);
					}
					else {
						Set<String> childrens = new HashSet<>();
						childrens.add(failingSetNode.toString());
						nodeChildrenMap.put(partialEmbedding.toString(), childrens);
					}
				}
				i++;
			}
		}
		return results;
	}
	
	private String checkChildNode(Map<String, Set<String>> failingSets, Map<String, Set<String>> nodeChildrenMap,
			Map<String, Integer> partialEmbedding, String nextExtendableNode) {
		// TODO Auto-generated method stub
		String child = "";
		for(String childNode : nodeChildrenMap.get(partialEmbedding.toString())) {
			if(!failingSets.get(childNode).contains(nextExtendableNode)) {
				child = childNode;
				break;
			}
		}
		return child;
	}

	private String convertEmbeddingToString(Map<String, Integer> sorted) {
		// TODO Auto-generated method stub
		StringBuilder sb = new StringBuilder();
		sb.append("S:8:");
		for(Entry<String, Integer> e : sorted.entrySet()) {
			sb.append(e.getKey().substring(1) + "," + e.getValue()+";");
		}
		sb.deleteCharAt(sb.length()-1);
		return sb.toString();
	}

	public void performDAFSubgraphMatching(String groundTruthFileName)
	{
		int queryNo = 0;
		try{
			
			//File groundTruthFolder = new File(groundTruthFolderPath);	
			File groundTruthFile = new File(groundTruthFolderPath+"/"+groundTruthFileName);
			BufferedReader bf_gtr = new BufferedReader(new FileReader(groundTruthFile));
			LinkedList<String> gtrResults = new LinkedList<>();
			String currLine = "";
			while ((currLine = bf_gtr.readLine()) != null)
            {
				if(!currLine.equals("") && !currLine.isEmpty())
					gtrResults.add(currLine);
            }
			bf_gtr.close();
			while(!gtrResults.isEmpty())
			{
				int TP = 0, FP = 0, FN = 0;
				HashSet<String> gtrSolution = new HashSet<>();
				String labelFileName = gtrResults.remove().split(":")[1].split("\\.")[0];
				String queryFileName = gtrResults.remove().split(":")[1];
				Integer results = Integer.parseInt(gtrResults.remove().split(":")[1]);
				System.out.println( "Query No: " + (++queryNo) );
				System.out.println("Graph File:"+labelFileName+".grf");
				System.out.println("Query Graph File: "+queryFileName);
				System.out.println("Results in Ground Truth: "+results);
				
				for(int i = 0;i<results;i++)
				{
					gtrSolution.add(gtrResults.remove());
				}
				
				long startTime = System.currentTimeMillis();//System.nanoTime();
				
				//naiveSubgraphMatch(queryFileName, labelFileName+".grf");
				
				Graph<Vertex, DefaultEdge> queryGraph = createQueryGraph(queryFileName);
				List<String> currentSolution = DAF(queryGraph, labelFileName);
				
				long queryTime = (System.currentTimeMillis()-startTime);
				System.out.println("Time taken in ms:"+queryTime);
				
				System.out.println("Ground Truth Solutions: "+gtrSolution.size());
				System.out.println("Query Solutions: "+currentSolution.size());
				
				for(String solution : currentSolution)
				{
					if(gtrSolution.contains(solution))
						TP++;
					else
						FP++;
				}
				FN = gtrSolution.size() - TP;

				currentSolution.clear();	
				
				System.out.println("True Positives:"+TP);
				System.out.println("False Positives:"+FP);
				System.out.println("False Negatives:"+FN);
				
				System.out.println("***************");
			}
		}
		catch (FileNotFoundException e) {
			// TODO: handle exception
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public void performDAFSubgraphMatchingByPruningFailingSets(String groundTruthFileName)
	{
		int queryNo = 0;
		try{
			
			//File groundTruthFolder = new File(groundTruthFolderPath);	
			File groundTruthFile = new File(groundTruthFolderPath+"/"+groundTruthFileName);
			BufferedReader bf_gtr = new BufferedReader(new FileReader(groundTruthFile));
			LinkedList<String> gtrResults = new LinkedList<>();
			String currLine = "";
			while ((currLine = bf_gtr.readLine()) != null)
            {
				if(!currLine.equals("") && !currLine.isEmpty())
					gtrResults.add(currLine);
            }
			bf_gtr.close();
			while(!gtrResults.isEmpty())
			{
				int TP = 0, FP = 0, FN = 0;
				HashSet<String> gtrSolution = new HashSet<>();
				String labelFileName = gtrResults.remove().split(":")[1].split("\\.")[0];
				String queryFileName = gtrResults.remove().split(":")[1];
				Integer results = Integer.parseInt(gtrResults.remove().split(":")[1]);
				System.out.println( "Query No: " + (++queryNo) );
				System.out.println("Graph File:"+labelFileName+".grf");
				System.out.println("Query Graph File: "+queryFileName);
				System.out.println("Results in Ground Truth: "+results);
				
				for(int i = 0;i<results;i++)
				{
					gtrSolution.add(gtrResults.remove());
				}
				
				long startTime = System.currentTimeMillis();//System.nanoTime();
				
				//naiveSubgraphMatch(queryFileName, labelFileName+".grf");
				
				Graph<Vertex, DefaultEdge> queryGraph = createQueryGraph(queryFileName);
				List<String> currentSolution = DAFFailingSets(queryGraph, labelFileName);
				
				long queryTime = (System.currentTimeMillis()-startTime);
				System.out.println("Time taken in ms:"+queryTime);
				
				System.out.println("Ground Truth Solutions: "+gtrSolution.size());
				System.out.println("Query Solutions: "+currentSolution.size());
				
				for(String solution : currentSolution)
				{
					if(gtrSolution.contains(solution))
						TP++;
					else
						FP++;
				}
				FN = gtrSolution.size() - TP;

				currentSolution.clear();	
				
				System.out.println("True Positives:"+TP);
				System.out.println("False Positives:"+FP);
				System.out.println("False Negatives:"+FN);
				
				System.out.println("***************");
			}
		}
		catch (FileNotFoundException e) {
			// TODO: handle exception
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public SimpleDirectedGraph<Vertex, DefaultEdge> createDAGEx5(){
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = new SimpleDirectedGraph<>(DefaultEdge.class);
		Vertex u1 = new Vertex("u1", "A");
		Vertex u2 = new Vertex("u2", "B");
		Vertex u3 = new Vertex("u3", "C");
		Vertex u4 = new Vertex("u4", "B");
		Vertex u5 = new Vertex("u5", "A");
		Vertex u6 = new Vertex("u6", "D");
		Vertex u7 = new Vertex("u7", "E");
		Vertex u8 = new Vertex("u8", "E");
		Vertex u9 = new Vertex("u9", "F");
		queryDAG.addVertex(u1);
		queryDAG.addVertex(u2);
		queryDAG.addVertex(u3);
		queryDAG.addVertex(u4);
		queryDAG.addVertex(u5);
		queryDAG.addVertex(u6);
		queryDAG.addVertex(u7);
		queryDAG.addVertex(u8);
		queryDAG.addVertex(u9);
		queryDAG.addEdge(u1, u2);
		queryDAG.addEdge(u1, u3);
		queryDAG.addEdge(u1, u4);
		queryDAG.addEdge(u2, u5);
		queryDAG.addEdge(u2, u6);
		queryDAG.addEdge(u3, u7);
		queryDAG.addEdge(u4, u8);
		queryDAG.addEdge(u6, u7);
		queryDAG.addEdge(u7, u9);
		return queryDAG;
	}
	
	private Map<Integer, Map<String, Set<Integer>>> buildCSAdjacencyListEx5(){
		Map<Integer, Map<String, Set<Integer>>> CSAdjacencyList = new HashMap<>();
		 Map<String, Set<Integer>> map1 = new HashMap();
		 map1.put("u1_u2", new HashSet<>(Arrays.asList(2, 3)));
		 map1.put("u1_u3", new HashSet<>(Arrays.asList(4, 5, 6)));
		 map1.put("u1_u4", new HashSet<>(Arrays.asList(2)));
		 CSAdjacencyList.put(1, map1);
		 
		 
		 Map<String, Set<Integer>> map2 = new HashMap();
		 map2.put("u2_u5", new HashSet<>(Arrays.asList(1, 7)));
		 map2.put("u2_u6", new HashSet<>(Arrays.asList(8)));
		 map2.put("u4_u8", new HashSet<>(Arrays.asList(12, 13, 14, 15)));
		 CSAdjacencyList.put(2, map2);
		 
		 Map<String, Set<Integer>> map3 = new HashMap();
		 map3.put("u2_u5", new HashSet<>(Arrays.asList(1, 7)));
		 map3.put("u2_u6", new HashSet<>(Arrays.asList(9)));
		 CSAdjacencyList.put(3, map3);
		 
		 Map<String, Set<Integer>> map4 = new HashMap();
		 map4.put("u3_u7", new HashSet<>(Arrays.asList(10)));
		 CSAdjacencyList.put(4, map4);
		 
		 Map<String, Set<Integer>> map5 = new HashMap();
		 map5.put("u3_u7", new HashSet<>(Arrays.asList(11)));
		 CSAdjacencyList.put(5, map5);
		 
		 Map<String, Set<Integer>> map6 = new HashMap();
		 map6.put("u3_u7", new HashSet<>(Arrays.asList(12)));
		 CSAdjacencyList.put(6, map6);
		 
		 Map<String, Set<Integer>> map8 = new HashMap();
		 map8.put("u6_u7", new HashSet<>(Arrays.asList(10 , 11)));
		 CSAdjacencyList.put(8, map8);
		 
		 Map<String, Set<Integer>> map9 = new HashMap();
		 map9.put("u6_u7", new HashSet<>(Arrays.asList(12)));
		 CSAdjacencyList.put(9, map9);
		 
		 Map<String, Set<Integer>> map10 = new HashMap();
		 map10.put("u7_u9", new HashSet<>(Arrays.asList(16, 17)));
		 CSAdjacencyList.put(10, map10);
		 
		 Map<String, Set<Integer>> map11 = new HashMap();
		 map11.put("u7_u9", new HashSet<>(Arrays.asList(16, 17, 18, 19, 20)));
		 CSAdjacencyList.put(11, map11);
		 
		 Map<String, Set<Integer>> map12 = new HashMap();
		 map12.put("u7_u9", new HashSet<>(Arrays.asList(19)));
		 CSAdjacencyList.put(12, map12);
		 
		return CSAdjacencyList;
	}
	
	public void testFailingSets(){
		
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = createDAGEx5();
		Map<String, Set<String>> ancestorsmap = getAncestors(queryDAG);
		System.out.println(ancestorsmap);
		Map<Integer, Map<String, Set<Integer>>> adjacencyList = buildCSAdjacencyListEx5();
		Map<String, Set<Integer>> CS = new HashMap<>();
		CS.put("u1", new HashSet<>(Arrays.asList(1)));
		List<Vertex> DAGorder = Arrays.asList(new Vertex("u1","A"), new Vertex("u2","B"), new Vertex("u6","D"), new Vertex("u5","A"), new Vertex("u3","C"), new Vertex("u7","E"), new Vertex("u9","F"), new Vertex("u4","B"), new Vertex("u8","E"));
		//List<String> results = backtrackFailingSet(new LinkedHashMap<>(), queryDAG, CS, DAGorder, new LinkedHashSet<>(), adjacencyList, ancestorsmap);
		List<String> results = backtrackFailingSets(new LinkedHashMap<>(), queryDAG, CS, DAGorder, new LinkedHashSet<>(), adjacencyList, ancestorsmap, new HashMap<>(), new HashMap<>());
		System.out.println(results);
	}
	
	public static void main(String[] args) throws IOException {
		initDB();
		DAFSubgraphMatchingNeo4J daf = new DAFSubgraphMatchingNeo4J();
		
		/*Graph<Vertex, DefaultEdge> queryGraph = daf.createQueryGraph("a_testquery.4.sub.grf");
		daf.DAF(queryGraph, "a_test");*/
		
		/*Graph<Vertex, DefaultEdge> queryGraph = daf.createQueryGraph("backbones_198L.8.sub.grf");
		daf.DAF(queryGraph, "ecoli_1ZTA");*/
		
		/*Graph<Vertex, DefaultEdge> queryGraph = daf.createQueryGraph("backbones_198L.8.sub.grf");
		daf.DAF(queryGraph, "human_2TGF");*/
		
		/*Graph<Vertex, DefaultEdge> queryGraph = daf.createQueryGraph("backbones_1EMA.8.sub.grf");
		daf.DAF(queryGraph, "backbones_1O54");*/
		
		//System.out.println(daf.createNeighboursMap("ecoli_1ZTA"));
		
		//daf.testFailingSets();
		
		/*Graph<Vertex, DefaultEdge> queryGraph = daf.createQueryGraph("backbones_198L.8.sub.grf");
		daf.DAFFailingSets(queryGraph, "ecoli_1ZTA");*/
		
		daf.performDAFSubgraphMatching("Proteins.8.gtr");
		
		//daf.performDAFSubgraphMatchingByPruningFailingSets("Proteins.8.gtr");
		
		shutdownDB();
	}
}
