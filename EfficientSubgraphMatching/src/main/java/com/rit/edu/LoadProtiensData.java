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

import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.unsafe.batchinsert.BatchInserter;
import org.neo4j.unsafe.batchinsert.BatchInserters;

public class LoadProtiensData {
	
	private static String DBPath = "ProtiensDB";
	private static String dataFolderPath = "Proteins/target";
	
	public static void loadProtienData()
	{
		BatchInserter inserter = null;
		try {
			inserter = BatchInserters.inserter(new File(DBPath));
			
			long startTime = System.currentTimeMillis();//System.nanoTime();
			File dataFolder = new File(dataFolderPath);
			String [] filesList = dataFolder.list();
			Arrays.sort(filesList);
			
			for(int fileIndex = 0;fileIndex<filesList.length; fileIndex++)
			{
				File dataFile = new File(dataFolderPath+"/"+filesList[fileIndex]);
				LinkedList<String> linesList = new LinkedList<>();
				BufferedReader bf = new BufferedReader(new FileReader(dataFile));
				String line = null;
				while ((line = bf.readLine()) != null)
	            {
					linesList.add(line);
	            }
				
				int nodeCounter = Integer.parseInt(linesList.getFirst());
				linesList.removeFirst();
				for(int i = 0;i < nodeCounter; i++)
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
				System.out.println("Data for file loaded :"+dataFile);
				bf.close();
				
			}
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
		//String query = "MATCH(N1:a_test) WHERE N1.id = '5' RETURN N1;";
		//String query = "MATCH(N1:a_test:C)-[:HASCONNECTION]-(N2) WHERE N1.id = 5 RETURN N2.id;";
		//String query = "RETURN EXISTS( (:a_test {id: '2'})-[:HASCONNECTION]-(:a_test {id: '9'}) )";
		String query = "MATCH(N1:backbones_1QMH) RETURN COUNT(N1) ;";
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

	public static void main(String[] args) {
		//loadProtienData();
		checkData();

	}

}
