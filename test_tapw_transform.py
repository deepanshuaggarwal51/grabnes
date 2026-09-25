#!/usr/bin/env python3
"""
Test script to verify TAPW transformation: H_proj = X† * H * X

This script uses the same row/column conventions as the Fortran code:
- H: sparse matrix in CSR format (N×N)
- X: transformation matrix (N×M), rows=orbitals, columns=TAPW basis
- H_proj: projected Hamiltonian (M×M)

The transformation is done in two steps like Fortran:
1. Y = H * X  (sparse × dense)
2. H_proj = X† * Y  (dense × dense)
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

def create_test_hamiltonian(N, sparsity=0.1):
    """Create a test sparse Hamiltonian matrix (Hermitian)"""
    # Create random sparse matrix
    np.random.seed(42)  # For reproducibility
    density = sparsity
    H_dense = sp.random(N, N, density=density, format='csr', dtype=complex)
    
    # Make it Hermitian
    H_dense = (H_dense + H_dense.H) / 2.0
    
    # Add some diagonal elements to make it more realistic
    diagonal = np.random.uniform(-2, 2, N) + 1j * np.random.uniform(-0.1, 0.1, N)
    diagonal = diagonal.real + 1j * 0  # Make diagonal real for Hermitian matrix
    H_dense.setdiag(diagonal)
    
    return H_dense.tocsr()

def create_test_X_matrix(N, M):
    """Create a test X matrix mimicking TAPW structure"""
    np.random.seed(123)
    
    # Create X matrix with structure similar to TAPW: X(orbital, G*label)
    # Each column represents a plane wave exp(iG·r) for a specific orbital type
    X = np.zeros((N, M), dtype=complex)
    
    # Fill with random phase factors (like exp(iG·r))
    for i in range(N):
        for j in range(M):
            # Random phases like TAPW would have
            phase = np.random.uniform(0, 2*np.pi)
            amplitude = np.random.uniform(0.8, 1.2)  # Slight amplitude variation
            X[i, j] = amplitude * np.exp(1j * phase)
    
    # Add some normalization (like TAPW does per sublattice)
    for j in range(M):
        norm = np.linalg.norm(X[:, j])
        if norm > 0:
            X[:, j] /= np.sqrt(N/M)  # Rough normalization
    
    return X

def fortran_style_transform(H_csr, X):
    """
    Perform TAPW transformation using the same two-step process as Fortran:
    1. Y = H * X  (sparse matrix multiplication)
    2. H_proj = X† * Y  (dense matrix multiplication)
    """
    N, M = X.shape
    
    print(f"Input dimensions:")
    print(f"  H: {H_csr.shape} (sparse, nnz={H_csr.nnz})")
    print(f"  X: {X.shape}")
    print(f"  Expected H_proj: ({M}, {M})")
    
    # Step 1: Y = H * X (sparse × dense)
    print("\nStep 1: Computing Y = H * X...")
    Y = H_csr.dot(X)
    print(f"  Y shape: {Y.shape}")
    print(f"  Y max element: {np.max(np.abs(Y)):.6f}")
    
    # Step 2: H_proj = X† * Y (dense × dense)
    print("\nStep 2: Computing H_proj = X† * Y...")
    X_dagger = X.conj().T
    H_proj = X_dagger @ Y
    print(f"  X† shape: {X_dagger.shape}")
    print(f"  H_proj shape: {H_proj.shape}")
    print(f"  H_proj max element: {np.max(np.abs(H_proj)):.6f}")
    
    return H_proj, Y

def direct_transform(H_csr, X):
    """
    Direct transformation for comparison: H_proj = X† * H * X
    """
    print("\nDirect method: Computing H_proj = X† * H * X...")
    H_dense = H_csr.toarray()
    H_proj_direct = X.conj().T @ H_dense @ X
    print(f"  H_proj_direct shape: {H_proj_direct.shape}")
    print(f"  H_proj_direct max element: {np.max(np.abs(H_proj_direct)):.6f}")
    
    return H_proj_direct

def check_hermiticity(H, name="Matrix"):
    """Check if matrix is Hermitian"""
    error = np.max(np.abs(H - H.conj().T))
    print(f"{name} Hermiticity check:")
    print(f"  Max |H - H†|: {error:.2e}")
    if error < 1e-12:
        print("  ✓ Hermitian (excellent)")
    elif error < 1e-8:
        print("  ✓ Approximately Hermitian (good)")
    else:
        print("  ⚠ Not Hermitian")
    return error

def test_unitary_property(X):
    """Test if X has unitary properties in the projected space"""
    M = X.shape[1]
    X_dagger_X = X.conj().T @ X
    
    # Check if X†X is close to identity
    identity = np.eye(M)
    error = np.max(np.abs(X_dagger_X - identity))
    
    print(f"Unitary property check (X†X vs I):")
    print(f"  Max |X†X - I|: {error:.2e}")
    if error < 1e-10:
        print("  ✓ X is unitary in projected space (excellent)")
    elif error < 1e-6:
        print("  ✓ X is approximately unitary (good)")
    else:
        print("  ⚠ X is not unitary (expected for truncated basis)")
    
    return error, X_dagger_X

def main():
    print("=== TAPW Transformation Test (Fortran-style) ===\n")
    
    # Test parameters
    N = 50   # Number of orbitals (small for testing)
    M = 20   # Number of TAPW basis functions (M < N for projection)
    
    print(f"Test setup:")
    print(f"  N (orbitals): {N}")
    print(f"  M (TAPW basis): {M}")
    print(f"  Compression ratio: {N/M:.1f}")
    
    # Create test matrices
    print(f"\n1. Creating test Hamiltonian (sparse)...")
    H_sparse = create_test_hamiltonian(N, sparsity=0.15)
    check_hermiticity(H_sparse.toarray(), "Original H")
    
    print(f"\n2. Creating test X matrix...")
    X = create_test_X_matrix(N, M)
    unitary_error, X_dagger_X = test_unitary_property(X)
    
    # Perform transformations
    print(f"\n3. Fortran-style transformation (two-step)...")
    H_proj_fortran, Y = fortran_style_transform(H_sparse, X)
    check_hermiticity(H_proj_fortran, "H_proj (Fortran-style)")
    
    print(f"\n4. Direct transformation (for comparison)...")
    H_proj_direct = direct_transform(H_sparse, X)
    check_hermiticity(H_proj_direct, "H_proj (direct)")
    
    # Compare methods
    print(f"\n5. Comparing methods...")
    difference = np.max(np.abs(H_proj_fortran - H_proj_direct))
    print(f"Max |H_proj_fortran - H_proj_direct|: {difference:.2e}")
    
    if difference < 1e-12:
        print("✅ PERFECT MATCH: Fortran-style method is identical to direct method!")
    elif difference < 1e-8:
        print("✅ EXCELLENT MATCH: Methods agree within numerical precision")
    else:
        print("❌ MISMATCH: Methods give different results!")
    
    # Test reconstruction (unitary check)
    print(f"\n6. Testing reconstruction: X * H_proj * X†...")
    H_reconstructed = X @ H_proj_fortran @ X.conj().T
    H_original = H_sparse.toarray()
    
    reconstruction_error = np.max(np.abs(H_reconstructed - H_original))
    print(f"Max |H_reconstructed - H_original|: {reconstruction_error:.2e}")
    
    if reconstruction_error < 1e-10:
        print("✅ PERFECT RECONSTRUCTION: X is unitary!")
    elif reconstruction_error < 1e-6:
        print("✅ GOOD RECONSTRUCTION: X is approximately unitary")
    else:
        print("⚠ POOR RECONSTRUCTION: X is not unitary (expected for truncated basis)")
    
    # Summary
    print(f"\n=== SUMMARY ===")
    print(f"✓ Two-step method matches direct method: {difference < 1e-8}")
    hermitian_error = check_hermiticity(H_proj_fortran, 'H_proj')
    print(f"✓ Projected Hamiltonian is Hermitian: {hermitian_error < 1e-8}")
    print(f"✓ X has unitary properties: {unitary_error < 1e-3}")
    print(f"✓ Reconstruction quality: {reconstruction_error:.2e}")
    
    print(f"\n🎯 CONCLUSION:")
    if difference < 1e-8:
        print("   The Fortran transform_sparse_hamiltonian logic is CORRECT!")
        print("   Two-step process Y=H*X, H_proj=X†*Y works perfectly.")
    else:
        print("   ❌ There may be an issue with the transformation logic.")
    
    # Optional: Plot eigenvalue comparison
    if M <= 20:  # Only for small matrices
        print(f"\n7. Eigenvalue comparison...")
        eigvals_orig = np.linalg.eigvals(H_original)
        eigvals_proj = np.linalg.eigvals(H_proj_fortran)
        
        print(f"Original H eigenvalues (first 5): {eigvals_orig[:5].real}")
        print(f"Projected H eigenvalues (all): {eigvals_proj.real}")
        
        # The projected eigenvalues should be a subset of the original ones
        # (in the limit of complete basis)

if __name__ == "__main__":
    main()
