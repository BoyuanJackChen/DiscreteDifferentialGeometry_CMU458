"use strict";

/**
 * @module Projects
 */
class ScalarPoissonProblem {
	/**
	 * This class solves a {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf scalar poisson problem} on a surface mesh.
	 * @constructor module:Projects.ScalarPoissonProblem
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} M The mass matrix of the input mesh.
	 * @property {number} totalArea The total surface area of the input mesh.
	 */
	constructor(geometry) {
		// index vertices
		this.vertexIndex = indexElements(geometry.mesh.vertices);

		// TODO: build the laplace and mass matrices, and compute total area
		this.A = geometry.laplaceMatrix(this.vertexIndex); 
		this.M = geometry.massMatrix(this.vertexIndex); 
		this.totalArea = geometry.totalArea(); 
	}

	/**
	 * Computes the solution of the poisson problem ΔΦ=ρ, or Ax = -M(rho - rhoBar), where A
	 * is the positive definite laplace matrix and M is the mass matrix.
	 * @method module:Projects.ScalarPoissonProblem#solve
	 * @param {module:LinearAlgebra.DenseMatrix} rho A scalar density of vertices of the input mesh.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	solve(rho) {    // rho for ρ. 
		var Mrho = this.M.timesDense(rho);
		var rhoBar = Mrho.sum() / this.totalArea;
		var rhoBarMatrix = DenseMatrix.constant(rhoBar, rho.nRows(), );

		var b = this.M.timesDense(rhoBarMatrix.minus(rho));
		let llt = this.A.chol();
		let x = llt.solvePositiveDefinite(b);

		return x; 
	}
}
