package com.rit.edu;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.jgrapht.graph.SimpleGraph;

public class DAFSubgraphMatching {
	
	public void DAF(Graph<Vertex, DefaultEdge> queryGraph, Graph<Vertex, DefaultEdge> dataGraph)
	{
		buildDAG(queryGraph, dataGraph);
	}
	
	public void buildDAG(Graph<Vertex, DefaultEdge> queryGraph, Graph<Vertex, DefaultEdge> dataGraph)
	{
		Vertex root = selectQueryRoot(queryGraph, dataGraph);
		List<List<Vertex>> verticesByLevel = traverseGraphinBFSOrder(root, queryGraph);
		SimpleDirectedGraph<Vertex, DefaultEdge> queryDAG = buildDAGFromVertices(verticesByLevel, queryGraph, dataGraph);
		
	}
	
	private SimpleDirectedGraph<Vertex, DefaultEdge> buildDAGFromVertices(List<List<Vertex>> verticesByLevel, Graph<Vertex, DefaultEdge> queryGraph, Graph<Vertex, DefaultEdge> dataGraph) {
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
			List<List<Vertex>> labelsGroup = groupAndSortVertices(lowerLevel, queryGraph, dataGraph);
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
		System.out.println("Query DAG" + queryDAG);
		return queryDAG;
	}

	/**
	 * Get the order of labels in the data graph (Most infrequent labels in graph come earlier)
	 * @param dataGraph
	 */
	private List<String> getLabelOrder(Graph<Vertex, DefaultEdge> dataGraph) {
		Map<String, Integer> labelCountsMap = new HashMap<>();
		for(Vertex v : dataGraph.vertexSet())
		{
			if(labelCountsMap.containsKey(v.label))
			{
				labelCountsMap.put(v.label, labelCountsMap.get(v.label) + 1);
			} else {
				labelCountsMap.put(v.label, 1);
			}
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
	 * Vertices by grouped by Labels. In each group the vertices are sorted in descending order of vertex degrees.  
	 * @param lowerLevel
	 * @param queryGraph
	 * @param dataGraph
	 */
	private List<List<Vertex>> groupAndSortVertices(List<Vertex> lowerLevel, Graph<Vertex, DefaultEdge> queryGraph, Graph<Vertex, DefaultEdge> dataGraph) {
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
		List<String> labelOrder = getLabelOrder(dataGraph);
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
	 * Sort vertices in descending order of their degree.
	 * @param verticesList
	 * @param queryGraph
	 * @return
	 */
	private List<Vertex> sortVerticesByDegree(List<Vertex> verticesList, Graph<Vertex, DefaultEdge> queryGraph) {
		verticesList.sort((Vertex v1, Vertex v2) -> queryGraph.edgesOf(v2).size() - queryGraph.edgesOf(v1).size());
		return verticesList;
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

	private Vertex selectQueryRoot(Graph<Vertex, DefaultEdge> queryGraph, Graph<Vertex, DefaultEdge> dataGraph) {
		Vertex root = null;
		float value = Float.MAX_VALUE;
		for(Vertex queryVertex : queryGraph.vertexSet())
		{
			float currValue = (float) getInitialCandidateSetOfVertex(queryVertex, queryGraph, dataGraph).size() / queryGraph.edgesOf(queryVertex).size();
			if (currValue < value) {
				root = queryVertex;
				value = currValue;
			}
		}
		return root;
	}
	
	public List<Vertex> getInitialCandidateSetOfVertex(Vertex queryVertex, Graph<Vertex, DefaultEdge> queryGraph, Graph<Vertex, DefaultEdge> dataGraph){
		List<Vertex> candidateSet = new ArrayList<Vertex>();
		for(Vertex dataVertex : dataGraph.vertexSet()) {
			if(queryVertex.label == dataVertex.label && (dataGraph.edgesOf(dataVertex).size() >= queryGraph.edgesOf(queryVertex).size())) {
				 candidateSet.add(dataVertex);
			}
		}
		return candidateSet;
	}

	public static Graph<Vertex, DefaultEdge> buildDataGraph()
	{
		Graph<Vertex, DefaultEdge> dataGraph = new SimpleGraph<>(DefaultEdge.class);
		Vertex v1 = new Vertex("v1", "A");
		Vertex v2 = new Vertex("v2", "A");
		Vertex v3 = new Vertex("v3", "B");
		Vertex v4 = new Vertex("v4", "B");
		Vertex v5 = new Vertex("v5", "C");
		Vertex v6 = new Vertex("v6", "C");
		Vertex v7 = new Vertex("v7", "C");
		Vertex v8 = new Vertex("v8", "B");
		Vertex v9 = new Vertex("v9", "C");
		Vertex v10 = new Vertex("v10", "D");
		Vertex v11 = new Vertex("v11", "D");
		Vertex v12 = new Vertex("v12", "B");
		dataGraph.addVertex(v1);
		dataGraph.addVertex(v2);
		dataGraph.addVertex(v3);
		dataGraph.addVertex(v4);
		dataGraph.addVertex(v5);
		dataGraph.addVertex(v6);
		dataGraph.addVertex(v7);
		dataGraph.addVertex(v8);
		dataGraph.addVertex(v9);
		dataGraph.addVertex(v10);
		dataGraph.addVertex(v11);
		dataGraph.addVertex(v12);
		
		dataGraph.addEdge(v1, v3);
		dataGraph.addEdge(v1, v4);
		dataGraph.addEdge(v1, v5);
		dataGraph.addEdge(v1, v6);
		dataGraph.addEdge(v1, v7);
		
		dataGraph.addEdge(v2, v8);
		dataGraph.addEdge(v2, v9);
		
		dataGraph.addEdge(v3, v5);
		dataGraph.addEdge(v3, v10);
		
		dataGraph.addEdge(v4, v5);
		dataGraph.addEdge(v4, v10);
		
		dataGraph.addEdge(v5, v10);
		
		dataGraph.addEdge(v6, v8);
		dataGraph.addEdge(v6, v10);
		
		dataGraph.addEdge(v7, v8);
		dataGraph.addEdge(v7, v10);
		
		dataGraph.addEdge(v8, v11);
		
		dataGraph.addEdge(v9, v11);
		dataGraph.addEdge(v9, v12);
		
		return dataGraph;
	}
	
	public static Graph<Vertex, DefaultEdge> buildQueryGraph()
	{
		Graph<Vertex, DefaultEdge> queryGraph = new SimpleGraph<>(DefaultEdge.class);
		
		Vertex u1 = new Vertex("u1", "A");
		Vertex u2 = new Vertex("u2", "B");
		Vertex u3 = new Vertex("u3", "C");
		Vertex u4 = new Vertex("u4", "D");
		queryGraph.addVertex(u1);
		queryGraph.addVertex(u2);
		queryGraph.addVertex(u3);
		queryGraph.addVertex(u4);
		
		queryGraph.addEdge(u1, u2);
		queryGraph.addEdge(u1, u3);
		
		queryGraph.addEdge(u2, u3);
		queryGraph.addEdge(u2, u4);
		
		queryGraph.addEdge(u3, u4);
	
		return queryGraph;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Graph<Vertex, DefaultEdge> dataGraph = buildDataGraph();
		
		Graph<Vertex, DefaultEdge> queryGraph = buildQueryGraph();
		System.out.println("Query Graph:" + queryGraph);
		new DAFSubgraphMatching().DAF(queryGraph, dataGraph);
		
 	}

}
