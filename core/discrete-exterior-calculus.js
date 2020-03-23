"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		// Takes vertices to dual faces
		let vertices = geometry.mesh.vertices;
		let numVertices = vertices.length;

		let T = new Triplet(numVertices, numVertices); // Like a double-key dictionary. 
		// Designed for accessing our matrix. 
		for (let v of vertices){
			let area = geometry.barycentricDualArea(v);
			let index = vertexIndex[v];

			T.addEntry(area, index, index);
		}
		return SparseMatrix.fromTriplet(T);
	} 
	// Each entry on the diagonal is the barycentric area. Vertices might be 1, 0.5 or 
	// sth else. Doesn't matter. 

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		// Takes vertices to dual faces
		let edges = geometry.mesh.edges;
		let numEdges = edges.length;

		let T = new Triplet(numEdges, numEdges); // Like a double-key dictionary. 
		// Designed for accessing our matrix. 
		for (let e of edges){

			let index = edgeIndex[e];
			let cotanB = geometry.cotan(e.halfedge);
			let cotanA = geometry.cotan(e.halfedge.twin);

			T.addEntry((cotanA+cotanB)/2, index, index);
		}
		return SparseMatrix.fromTriplet(T);
	}
	// Key: we can compute the ratio of the dual edge and original edge
	// with the cotangent of two angles! See why later

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		let faces = geometry.mesh.faces;
		let numFaces = faces.length;
		let T = new Triplet(numFaces, numFaces);

		for (let f of faces) {
			let index = faceIndex[f];
			let area = geometry.area(f);
			T.addEntry(1/area, index, index);
		}
		return SparseMatrix.fromTriplet(T);
	}
	// Each entry is the 1/Area of the primal. We want to get areas into 1. 
	// ??? Why do we need this??????????

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		let edges = geometry.mesh.edges;
		let vertices = geometry.mesh.vertices; 
		// The number of cols is the number of vertices; the number of rows
		// is the number of edges. Parameter: rows → collumns
		let T = new Triplet(edges.length, vertices.length); // |E| x |V|

		// Let tail be the starting point; head is the head of the arrow
		for (let e of edges) {
			let rowIndex = edgeIndex[e];
			let headCol = vertexIndex[e.halfedge.vertex];
			let tailCol = vertexIndex[e.halfedge.twin.vertex];
			T.addEntry(-1, rowIndex, tailCol);
			T.addEntry(1, rowIndex, headCol);
		}
		// For each row, plug in -1 for head; 1 for tail

		return SparseMatrix.fromTriplet(T); 
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		let edges = geometry.mesh.edges;
		let faces = geometry.mesh.faces; 

		let T = new Triplet(faces.length, edges.length); // |F| x |E|

		for (let f of faces) {
			let rowIndex = faceIndex[f];
			for (let h of f.adjacentHalfedges()){
				let colIndex = edgeIndex[h.edge];

				let entry = 1; // Default: the edge and halfedge agree in orientation
				if (h != h.edge.halfedge) {
					entry = -1;
				}
				T.addEntry(entry, rowIndex, colIndex);
			}
		}

		return SparseMatrix.fromTriplet(T); 
	}
}
