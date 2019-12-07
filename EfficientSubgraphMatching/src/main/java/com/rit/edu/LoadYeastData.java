package com.rit.edu;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.unsafe.batchinsert.BatchInserter;
import org.neo4j.unsafe.batchinsert.BatchInserters;

public class LoadYeastData {

	private static String DBPath = "YeastDB";
	private static String yeastDataFile = "HumanAndYeast/yeast.igraph";
	private static String queryDataFile = "HumanAndYeast/yeast.igraph";
	
	public static void loadYeastData()
	{
		BatchInserter inserter = null;
		try {
			inserter = BatchInserters.inserter(new File(DBPath));
			
			long startTime = System.currentTimeMillis();//System.nanoTime();
			
			File dataFile = new File(yeastDataFile);
			LinkedList<String> linesList = new LinkedList<>();
			BufferedReader bf = new BufferedReader(new FileReader(dataFile));
			String line = null;
			while ((line = bf.readLine()) != null)
            {
				linesList.add(line);
            }
			
			for(String l : linesList) {
				if(l.contains("t #"))
					continue;
				if(l.startsWith("v")) {
					String [] arr = l.split(" ");
					//System.out.println("Vertex " + arr[1]);
					long id = (long) Integer.parseInt( arr[1]);
					//String labels [] = Arrays.copyOfRange(arr, 2, arr.length);
					Label [] labels = new Label[arr.length - 2];
					for(int i = 0; i< labels.length; i++) {
						labels[i] = Label.label(arr[i + 2]);
					}
					Map<String,Object> attributes = new HashMap<String, Object>();
					attributes.put("id", arr[1]);
					inserter.createNode(id, attributes, labels);
					//System.out.println(Arrays.toString(labels));
					
				}
				if(l.startsWith("e")) {
					String [] relationArr = l.split(" ");
					long id1 = (long) Integer.parseInt(relationArr[1]);
					long id2 = (long) Integer.parseInt(relationArr[2]);
					inserter.createRelationship(id1, id2, RelationshipType.withName("HASCONNECTION"), null);
					inserter.createRelationship(id2, id1, RelationshipType.withName("HASCONNECTION"), null);
				}
			}
			/*for(int i = 0;i < nodeCounter; i++)
			{
				String nodeArr [] = linesList.getFirst().split(" ");
				long id = (long) Integer.parseInt(nodeArr[0]) + fileIndex * 100000 ;
				String label1 = nodeArr[1];
				String label2 = filesList[fileIndex].split(".grf")[0];
				Map<String,Object> attributes = new HashMap<String, Object>();
				attributes.put("id", nodeArr[0]);
				inserter.createNode(id, attributes, Label.label(label1), Label.label(label2));
				linesList.removeFirst();
			}
			while(!linesList.isEmpty())
			{
				int relationCounter = Integer.parseInt(linesList.getFirst());
				linesList.removeFirst();
				for(int i = 0; i< relationCounter;i++)
				{
					String relationArr [] = linesList.getFirst().split(" ");
					linesList.removeFirst();
					long id1 = (long) Integer.parseInt(relationArr[0]) + fileIndex * 100000 ;
					long id2 = (long) Integer.parseInt(relationArr[1]) + fileIndex * 100000 ;
					inserter.createRelationship(id1, id2, RelationshipType.withName("HASCONNECTION"), null);
				}
			}
			System.out.println("Data for file loaded :"+dataFile);*/
			
			bf.close();
			long queryTime = (System.currentTimeMillis()-startTime);
			System.out.println("Total Time in ms :"+queryTime);
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		finally {
			inserter.shutdown();
			System.out.println("Database shutdown complete.");
		}
	}
	
	private static void checkData()
	{
		//String query = "MATCH(N1:a_test) RETURN COUNT(N1) ;";
		//String query = "MATCH(N1:a_test) RETURN N1.id ;";
		String query = "MATCH(N1) WHERE N1.id = '2' RETURN labels(N1);";
		//String query = "MATCH(N1:a_test:C)-[:HASCONNECTION]-(N2) WHERE N1.id = 5 RETURN N2.id;";
		//String query = "RETURN EXISTS( (:a_test {id: '2'})-[:HASCONNECTION]-(:a_test {id: '9'}) )";
		GraphDatabaseFactory dbFactory = new GraphDatabaseFactory();
		GraphDatabaseService dbService = dbFactory.newEmbeddedDatabase(new File(DBPath));
		Result rs = dbService.execute(query);
		while(rs.hasNext())
		{
			Map<String,Object> next = rs.next();
			for(Entry<String, Object> entry : next.entrySet())
			{
				System.out.println(entry.getKey() +" : "+entry.getValue().toString());
			}
		}
		rs.close();
		dbService.shutdown();
	}
	
	
	public Graph<Vertex, DefaultEdge> createQueryGraph(int queryNo) throws IOException {
		//SimpleGraph<Vertex, DefaultEdge> queryGraph = new SimpleGraph<Vertex, DefaultEdge>(DefaultEdge.class);
		Graph<Vertex, DefaultEdge> queryGraph = new SimpleGraph<>(DefaultEdge.class);
		File queryFile = new File(queryDataFile);
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

	public static void main(String[] args) {
		//loadYeastData();
		checkData();

	}

}
