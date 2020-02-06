"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

    /** This class implements various operators (e.g. boundary, star, link) on a mesh.
     * @constructor module:Projects.SimplicialComplexOperators
     * @param {module:Core.Mesh} mesh The input mesh this class acts on.
     * @property {module:Core.Mesh} mesh The input mesh this class acts on.
     * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
     * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
     */
    constructor(mesh) {
        this.mesh = mesh;
        this.assignElementIndices(this.mesh);

        this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
        this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
    }

    /** Assigns indices to the input mesh's vertices, edges, and faces
     * @method module:Projects.SimplicialComplexOperators#assignElementIndices
     * @param {module:Core.Mesh} mesh The input mesh which we index.
     */
    assignElementIndices(mesh) {
		this.mesh.indexElements();
    }

    /** Returns the vertex-edge adjacency matrix of the given mesh.
     * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
     * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
     * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
     */
    buildVertexEdgeAdjacencyMatrix(mesh) {
		// Initialize the vertex-edge matrix. First layer is edge; second layer is vertices.
		var matrix = new Array(mesh.edges.length)
		for (var i = 0; i < mesh.edges.length; i++) {
			matrix[i] = new Array(mesh.vertices.length).fill(0);
		}
		// For each vertex, iterate through the halfedges. 
        for (var i = 0; i < mesh.vertices.length; i++) {
            var vhi = new VertexHalfedgeIterator(mesh.vertices[i].halfedge);
			for(let thisHalfEdge of vhi){
				matrix[thisHalfEdge.edge.index][mesh.vertices[i].index] = 1
			}
        }
		return matrix;
    }

    /** Returns the edge-face adjacency matrix.
     * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
     * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
     * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
     */
    buildEdgeFaceAdjacencyMatrix(mesh) {
		// Initialize the vertex-edge matrix. First layer is faces; second layer is edges.
		var matrix = new Array(mesh.faces.length)
		for (var i = 0; i < mesh.faces.length; i++) {
			matrix[i] = new Array(mesh.edges.length).fill(0);
		}
		// For each face, iterate through the halfedges. 
        for (var i = 0; i < mesh.faces.length; i++) {
            var fhi = new FaceHalfedgeIterator(mesh.faces[i].halfedge);
			for(let thisHalfEdge of fhi){
				matrix[mesh.faces[i].index][thisHalfEdge.edge.index] = 1
			}
        }
		return matrix;
    }

    /** Returns a column vector representing the vertices of the
     * given subset.
     * @method module:Projects.SimplicialComplexOperators#buildVertexVector
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
     *  vertex i is in the given subset and 0 otherwise
     */
    buildVertexVector(subset) {
		var vector = new Array(mesh.vertices.length).fill(0);
		console.log('Subset number of vertices: ');
		console.log(subset.vertices.size);
		for(let v_index of subset.vertices) {
			// console.log(Object.prototype.toString.call(v_index));
			vector[v_index] = 1;
		}
		return vector;
    }

    /** Returns a column vector representing the edges of the
     * given subset.
     * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
     *  edge i is in the given subset and 0 otherwise
     */
    buildEdgeVector(subset) {
		var vector = new Array(mesh.edges.length).fill(0);
		for(let e_index of subset.edges) {
			vector[e_index] = 1;
		}
		return vector;
    }

    /** Returns a column vector representing the faces of the
     * given subset.
     * @method module:Projects.SimplicialComplexOperators#buildFaceVector
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
     *  face i is in the given subset and 0 otherwise
     */
    buildFaceVector(subset) {
		var vector = new Array(mesh.faces.length).fill(0);
		for(let f_index of subset.faces) {
			vector[f_index] = 1;
		}
		return vector;
    }

	// Debugging function
	print2DMatrix(matrix) {
  		for(var z = 0; z < matrix.length; z++) {
			for(var i = 0; i < matrix[0].length; i++) {
    			console.log(matrix[z][i]);
  			}
		}
	}

	// Debugging function
	print1Entries(matrix) {
		var accum = 0;
		console.log("The 1 entries are: ");
		for(var z = 0; z < matrix.length; z++) {
			for(var i = 0; i < matrix[0].length; i++) {
				if (matrix[z][i]==1) {
    				console.log(z,i);
					accum++;
				}
  			}
		}
		console.log(accum/2);
		console.log("End of 1 entries");
	}

    /** Returns the star of a subset.
     * @method module:Projects.SimplicialComplexOperators#star
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:Core.MeshSubset} The star of the given subset.
     */
    star(subset) {
		// The star should at least have all the original components
		var result = MeshSubset.deepCopy(subset);  
		var vertexVector = this.buildVertexVector(subset);
		var edgeVector = this.buildEdgeVector(subset);

		// Process vertices
		for (var i = 0; i < vertexVector.length; i++) {
			if (vertexVector[i]==1) {
				// Add related edges
				var vei = new VertexEdgeIterator(this.mesh.vertices[i].halfedge);
				for (let edge of vei) {
					result.addEdge(edge.index);
				}
				// Add related faces
				var vfi = new VertexFaceIterator(this.mesh.vertices[i].halfedge);
				for (let face of vfi) {
					result.addFace(face.index);
				}
			}
		}
		// Process edges
		for (var i = 0; i < edgeVector.length; i++) {
			if (edgeVector[i]==1) {
				// Add related faces
				for (var j = 0; j < this.A1.length; j++) {
					if (this.A1[j][i]==1){
						result.addFace(this.mesh.faces[j].index);
					}
				}
			}
		}
		// Finally, Faces don't generate new components for star.
		// So we don't need to process faces. 
        return result;
    }

    /** Returns the closure of a subset.
     * @method module:Projects.SimplicialComplexOperators#closure
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:Core.MeshSubset} The closure of the given subset.
     */
    closure(subset) {
		// The closure should at least have all the original components
		var result = MeshSubset.deepCopy(subset); 
		var edgeVector = this.buildEdgeVector(subset); 
		var faceVector = this.buildFaceVector(subset);
		// Vertices do not lead to extra elements in closre. 

        // Process edges by adding two related vertices
		for (var i = 0; i < edgeVector.length; i++) {
			if (edgeVector[i]==1) {
				var v1_index = this.mesh.edges[i].halfedge.vertex.index;
				result.addVertex(v1_index);
				var v2_index = this.mesh.edges[i].halfedge.next.vertex.index;
				result.addVertex(v2_index); 
			}
		}
		// Process faces by adding related edges and vertices
		for (var i = 0; i < faceVector.length; i++) {
			if (faceVector[i]==1) {
				var vertices  = this.mesh.faces[i].adjacentVertices();
				var edges = this.mesh.faces[i].adjacentEdges();
				for (let v of vertices) {
					result.addVertex(v.index);
				}
				for (let e of edges) {
					result.addEdge(e.index);
				}
			}
		}
        return result; 
    }

    /** Returns the link of a subset.
     * @method module:Projects.SimplicialComplexOperators#link
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:Core.MeshSubset} The link of the given subset.
     */
    link(subset) {
		var left  = this.closure(this.star(subset));
		var right = this.star(this.closure(subset));
		left.deleteVertices(right.vertices);
		left.deleteEdges(right.edges);
		left.deleteFaces(right.faces);
        return left; 
    }

    /** Returns true if the given subset is a subcomplex and false otherwise.
     * @method module:Projects.SimplicialComplexOperators#isComplex
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
     */
    isComplex(subset) {
		return (this.closure(subset).equals(subset));
    }

    /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
     * @method module:Projects.SimplicialComplexOperators#isPureComplex
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
     */
    isPureComplex(subset) {
		if (this.isComplex(subset)) {
			if (subset.faces.size > 0)
				return 2;
			else if (subset.edges.size > 0)
				return 1;
			else if (subset.vertices.size > 0)
				return 0;
			else 
				return -1;
		} 
		else 
			return -1;
    }

    /** Returns the boundary of a subset.
     * @method module:Projects.SimplicialComplexOperators#boundary
     * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
     * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
     */
    boundary(subset) {
		var result = new MeshSubset();
		// Add the border edges by going through faces
		for (let f_index of subset.faces) {
			var fei = this.mesh.faces[f_index].adjacentEdges();
			for (let edge of fei) {
				if (edge.onBoundary()) {
					result.addEdge(edge.index);
				}
				else if (!(subset.faces.has(edge.halfedge.face) 
					&& subset.faces.has(edge.halfedge.twin.face))) {
					result.addEdge(edge.index);
				}
			}
		}
		// Then add the vertices attached to these edges
		for (let e_index of result.edges) {
			var v_index1 = this.mesh.edges[e_index].halfedge.vertex.index;
			var v_index2 = this.mesh.edges[e_index].halfedge.twin.vertex.index;
			result.addVertex(v_index1);
			result.addVertex(v_index2);
		}
        return result; 
    }
}
