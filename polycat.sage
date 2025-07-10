# ------------------------------------------------------------------------------
#
# Copyright (C) 2025 Arnau Padrol <arnau.padrol@ub.edu>
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
#          https://www.gnu.org/licenses/
#
# ------------------------------------------------------------------------------
#
# This module requires SageMath (https://www.sagemath.org/)
# for matrix, vector, polyhedron, and graph operations.
# It is not compatible with pure Python or standard scientific libraries alone.
#
# ------------------------------------------------------------------------------


def ith_projection(frame, i, project_on_span=False):
    r"""
    Gives the matrix of the ith projection associated to a frame.
   
    INPUT:
    
    - ``frame``: a list of d linearly independent vectors, each with d coordinates (i.e., a basis)
    - ``i``: the dimension onto which to project (i <= d)
    - ``project_on_span``: bool (default False)
       If False, returns the i x d projection matrix onto the span of the first i vectors, expressed in frame coordinates.
       If True, returns the projection on the span of the frame corresponding to the first i vectors, in ambient coordinates.


    OUTPUT:
    
    - The requested projection matrix.
      
    EXAMPLE:
    
        sage: frame = [vector([1,0,0]), vector([0,1,0]), vector([0,0,1])]
        sage: ith_projection(frame, 2)
        [1 0 0]
        [0 1 0]
    """
    

    d = len(frame)
    
    # Check i is in the valid range
    if not (1 <= i <= d):
        raise ValueError(f"Parameter i must satisfy 1 <= i <= {d}, got i={i}")
    
    # Check all vectors have length d
    for idx, v in enumerate(frame):
        if len(v) != d:
            raise ValueError(f"Vector at index {idx} has length {len(v)}, expected {d}")
    
    M = Matrix(frame).transpose()
    
    # Check M is invertible (i.e., frame vectors are linearly independent)
    if M.det() == 0:
        raise ValueError("Frame vectors are not linearly independent")
    

    if not project_on_span:
        # Projection in R^i
        I = block_matrix([[identity_matrix(i) , zero_matrix(i, d-i)]])
        Proj=I*M.inverse()
    else:
        # Projection on the span of the frame
        I = block_matrix([
            [identity_matrix(i), zero_matrix(i, d - i)],
            [zero_matrix(d - i, i), zero_matrix(d - i)]
        ])
        Proj = M * I * M.inverse()
    
    return Proj



def system_of_projections(frame, project_on_span=False):
    """
    Computes the full system of projection matrices associated to a frame.
    
    INPUT:
        - ``frame``: a list of d linearly independent vectors forming a basis.
        - ``project_on_span``: bool (default False)
       If False, returns the i x d projection matrix onto the span of the first i vectors, expressed in frame coordinates.
       If True, returns the projection on the span of the frame corresponding to the first i vectors, in ambient coordinates.

    
    OUTPUT:
        - A list S of length d+1 where S[i] is the projection matrix onto the span of
          the first i vectors of the frame (with S[0] being the zero projection), in frame coordinates or ambient coordinates
    
    NOTES:
        - Uses the ith_projection function to compute each projection matrix.
    """
    d = len(frame)
    S = []
    for i in range(d + 1):
        S.append(ith_projection(frame, i, project_on_span=project_on_span))
        
    return S



def canonical_frame(d):
    """
    Returns the canonical frame of dimension d.
    
    INPUT:
        - ``d``: positive integer, the dimension.
        
    OUTPUT:
        - A list of d vectors forming the standard frame of R^d.
    """
    return identity_matrix(d).columns()


def canonical_system_of_projections(d):
    """
    Returns the system of projections associated with the canonical frame of R^d.
    
    INPUT:
        - ``d``: positive integer, the dimension.
        
    OUTPUT:
        - A list of d+1 projection matrices corresponding to the canonical frame.
    """
    return system_of_projections(canonical_frame(d))



def admissible(frame,Pol,verbose=False):
    
    for i in range(1, Pol.dim()):
        for F in Pol.faces(i):
            F=F.as_polyhedron()
            Proj=ith_projection(frame,i)
            if (Proj*F).dim()!= i:
                if verbose: 
                    print('the basis is not admissible because of this face', F, F.vertices(), ' projects to ', Proj*F, (Proj*F).vertices(),)
                return False
    return True

def is_admissible(frame, Pol,verbose=False):
    """
    Checks if a given frame is admissible with respect to a polytope Pol.

    INPUT:
        - ``frame``: a list of d linearly independent vectors forming a basis.
        - ``Pol``: a polytope in R^d

    OUTPUT:
        - Returns True if frame is admissible for Pol, False otherwise, printing details of the first violating face.

    """
    

    d = len(frame)
    
    # Check all vectors have length d
    for idx, v in enumerate(frame):
        if len(v) != d:
            raise ValueError(f"Vector at index {idx} has length {len(v)}, expected {d}")
    
    M = Matrix(frame).transpose()
    
    # Check M is invertible (i.e., frame vectors are linearly independent)
    if M.det() == 0:
        raise ValueError("Frame vectors are not linearly independent")
        
    # Check P is in R^d  
    if P.ambient_dim()!=d:
        raise ValueError(f"The polytope's ambient dimension is {P.ambient_dim()}, expected {d}")
    
    for i in range(1, Pol.dim()):
        Proj = ith_projection(frame, i)
        for F in Pol.faces(i):
            F = F.as_polyhedron() 
            projected_F = Proj * F
            
            if projected_F.dim() != i:
                if verbose:
                    print("The basis is not admissible because of this face:", F)
                    print("Face vertices:", F.vertices())
                    print("Projects to:", projected_F)
                    print("Projected vertices:", projected_F.vertices())
                return False
    
    return True


def source_target(Pol, i, frame='canonical'):
    """
    Returns the i-dimensional source and target of a polytope Pol
    with respect to a given frame (default is the canonical frame).

    INPUT:
        - Pol: a polytope 
        - i: integer dimension 
        - frame: a list of vectors forming a basis, or 'canonical' to use the standard frame.

    OUTPUT:
        - A tuple (source, target), where each is a set of polytopes corresponding to the
          source and target i-dimensional faces of Pol.

    """
    if frame == 'canonical':
        frame = canonical_frame(Pol.ambient_dimension())
        
    d = len(frame)
    if not (0 <= i < d):
        raise ValueError(f"Parameter i must satisfy 0 <= i < {d}, got i={i}")
    
    Proj = ith_projection(frame, i + 1)
    ProjPol = Proj * Pol
    orientation = Proj * frame[i]
    
#    print(orientation)
    
    target = set()
    source = set()
    
    for F in ProjPol.facets():
        nv = F.normal_cone().rays()[0].vector()
        if nv * orientation > 0:
            target.add(F.as_polyhedron())
        elif nv * orientation < 0:
            source.add(F.as_polyhedron())
    
    return (source, target)


def source(Pol, i, frame='canonical'):
    """
    Returns the i-dimensional source faces of Pol for the given frame.

    OUTPUT:
        - A set of source polyhedra (faces).
    """
    return source_target(Pol, i, frame)[0]


def target(Pol, i, frame='canonical'):
    """
    Returns the i-dimensional target faces of Pol for the given frame.

    OUTPUT:
        - A set of target polyhedra (faces).
    """
    return source_target(Pol, i, frame)[1]



def cell_k_string(Pol, k, frame='canonical', pure=False, as_labels=True):
    """
    Computes the directed graph containing the cellular k-string of a polytope Pol.

    Vertices of the graph are faces of Pol, and there is a directed edge F -> G if
    the target of F intersects the source of G.

    INPUT:
        - Pol: a polytope object with methods `.ambient_dimension()` and `.faces(dim)`.
        - k: integer dimension for the source/target computation.
        - frame: frame used for source_target, default 'canonical'.
        - pure: boolean, if True only consider faces of dimension k+1,
                if False consider faces of dimension from k+1 up to ambient dimension.
        - as_labels: boolean, if True vertices are labeled by sets of vertex indices,
                     otherwise by the face objects themselves.

    OUTPUT:
        - A directed graph (DiGraph) whose vertices are faces or their labels,
          and edges correspond to adjacency via source-target intersections.
    """
    d = Pol.ambient_dimension()
    if not (0 <= k < d):
        raise ValueError(f"Parameter k must satisfy 0 <= k < {d}, got k={k}")

    facelist = []
    
    # Select dimensions to iterate over
    dims_to_consider = [k + 1] if pure else list(range(k + 1, d))

              
    for dim in dims_to_consider:
        for F in Pol.faces(dim):
            source, target = source_target(F.as_polyhedron(), k, frame)
            facelist.append((F, source, target))

            
    # Precompute labels if needed to avoid recomputing inside nested loops
    if as_labels:
        labeled_faces = [(Set(F.ambient_V_indices()), source, target) for F, source, target in facelist]
    else:
        labeled_faces = facelist
        
        
    D = DiGraph()
    
    # Add vertices
    for Fi, _, _ in labeled_faces:
        D.add_vertex(Fi)
        
    # Add edges according to intersection of source/target sets
    for i, (Fi, source_i, target_i) in enumerate(labeled_faces):
        for j in range(i + 1, len(labeled_faces)):
            Fj, source_j, target_j = labeled_faces[j]

            if target_i.intersection(source_j):
                D.add_edge(Fi, Fj)
            if target_j.intersection(source_i):
                D.add_edge(Fj, Fi)

    return D

    