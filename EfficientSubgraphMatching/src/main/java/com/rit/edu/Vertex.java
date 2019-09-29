package com.rit.edu;

public class Vertex {
	
	String id;
	String label;
	
	public Vertex(String id, String label) {
		super();
		this.id = id;
		this.label = label;
	}
	
	@Override
	public String toString() {
		return id + ":" + label;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Vertex) {
			Vertex v = (Vertex) obj;
			if( this.id.equals(v.id) && this.label.equals(label)) {
				return true;
			}
		}
		return false;
	}
}
