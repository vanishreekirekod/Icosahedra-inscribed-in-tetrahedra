with(SimplicialSurfaceEmbeddings):
with(combinat):
#with(LinearAlgebra[Modular]):
with(ListTools):
with(LinearAlgebra):


#converts list of lists to set of sets.
convertLLtoSS := proc(S) 
local i, j, L; L := {}; 
for i to nops(S) do 
    L := L union {convert(S[i], set)}; 
end do; 
return L; 
end proc:

#converts set of sets to list of lists
convertSStoLL := proc(S) 
local i, j, L; L := []; 
for i to nops(S) do 
    L := [op(L), convert(S[i], list)]; 
end do; 
return L; 
end proc:

#tetrahedra
tetrahedron:=proc(var::list)
local surf,coords,ev,vert,ori,tet_ed,numV,s1,s2,co;
surf:=NewSurface():
AddFace(surf,[1,2,3],"lengths"=[var[3],var[2],var[1]]):  
AttachFace(surf,[1,2],4,"lengths"=[var[3],var[1],var[2]]):  
coords:=CoordinateMatrix(surf, 1,"listlist"=true):  
SquareDistance(coords[3],coords[4]): 
ev:=evala(`%`): 
var[3]^2-ev:
s1,s2:=solve(`%`,_t[4]):  
co:=subs(_t[4]=s1,coords):
numV:=4:
vert:=[1,2,3,4]:
ori:=[[1, 2, 3], [2, 1, 4], [3, 2, 4], [4, 1, 3]]:
tet_ed:=[[1, 2], [2,3], [1,3], [1,4], [2,4], [3,4]]:
return simplify(co),ori,vert;
end proc:


#octahedra from tetrahedron
octahedron_from_tetrahedron:=proc(var::list)
local co_O, ii, fac_O,vert_O,surf,coords,ev,s1,s2,co,numV,ori,tet_ed,vert;
surf:=NewSurface():
AddFace(surf,[1,2,3],"lengths"=[var[3],var[2],var[1]]):  
AttachFace(surf,[1,2],4,"lengths"=[var[3],var[1],var[2]]):  
coords:=CoordinateMatrix(surf, 1,"listlist"=true):  
SquareDistance(coords[3],coords[4]): 
ev:=evala(`%`): 
var[3]^2-ev:
s1,s2:=solve(`%`,_t[4]):  
co:=subs(_t[4]=s1,coords):
numV:=4:
vert:=[1,2,3,4]:
ori:=[[1, 2, 3], [2, 1, 4], [3, 2, 4], [4, 1, 3]]:
tet_ed:=[[1, 2], [2,3], [1,3], [1,4], [2,4], [3,4]]:
co_O:=[]:
for ii from 1 to nops(tet_ed) do  
    co_O:=[op(co_O),(co[tet_ed[ii][2]]+co[tet_ed[ii][1]])/(2)]; 
end do:
fac_O:=[[1,3,4],[5,2,1],[3,2,6],[4,6,5],[1,2,3],[4,5,1],[6,2,5],[4,3,6]], [1, 2, 3, 4, {1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {1, 3, 4}]:
vert_O:=[1,2,3,4,5,6]:
return simplify(co_O),fac_O[1],vert_O;
end proc:



goldenratio:=proc(c1,c2,dist,var_ratio::list)
local i,j,C;
C:=[((var_ratio[2]*c1[1])+(var_ratio[1]*c2[1]))/(var_ratio[1]+var_ratio[2]), ((var_ratio[2]*c1[2])+(var_ratio[1]*c2[2]))/(var_ratio[1]+var_ratio[2]), ((var_ratio[2]*c1[3])+(var_ratio[1]*c2[3]))/(var_ratio[1]+var_ratio[2])]; 
return C;
end proc:


icosahedron:=proc(co_O,fac_O,vert_O,var::list)
local i,j,oct_ed,oct_edel,co_I,temp_ed,vert_I,fac_I,temp_oct_ed,ed1,temp,coo,m,n; 
oct_ed:=[[2,1],   [1,3], [4,1], [3,2], [5,2], [1,5], [3,4], [6,3], [2,6], [6,5], [5,4], [4,6]];
oct_edel:=[var[1]^2/4, var[2]^2/4, var[1]^2/4, var[3]^2/4, var[3]^2/4, var[2]^2/4, var[3]^2/4, var[2]^2/4, var[1]^2/4, var[2]^2/4, var[3]^2/4, var[1]^2/4];
fac_I:=[[12, 10, 11], [11, 10, 6], [6, 10, 5], [5, 1, 6], [1, 3, 6], [11, 6, 3], [3, 7, 11], [3, 2, 7], [1, 2, 3], [2, 1, 4], [1, 5, 4], [4, 5, 9], [9, 5, 10], [9, 10, 12], [12, 8, 9], [4, 9, 8], [2, 4, 8], [7, 2, 8], [7, 8, 12], [7, 12, 11]];
co_I:=[];
vert_I:=[seq(i,i=1..12)];
temp_oct_ed:={};  
for i from 1 to nops(oct_ed) do 
    coo:=simplify(goldenratio(co_O[oct_ed[i][1]],co_O[oct_ed[i][2]],oct_edel[i],[var[4],var[5]])); 
    co_I:=[op(co_I),coo]; 
end do;
return simplify(co_I),fac_I,vert_I;
end proc:


d:=tetrahedron([a,b,c]):
dd:=octahedron_from_tetrahedron([a,b,c]):
ddd:=icosahedron(dd[1],dd[2],dd[3],[a,b,c,m,n]): ddd[2],ddd[3];


####Automorphism Group
InvPer := proc(a::list)
local j,ai;
ai := map(i->0,a):
for j from 1 to nops(a) do
  ai[a[j]] := j;
od;
ai;
end:

ProPer := proc(a::list, b::list)
local c,j;
c := map(x->0, [$1..nops(a)]):
for j from 1 to nops(a) do
  c[j] := a[b[j]];
od;
c;
end:

Bahn := proc(G::list, m, omega)
local B, A,aa;
B := {};
A := {m};
while A <> {} do
  aa := A[1]; #print(op(map(g->omega(g,aa),G)),"next"):
  A := {op(A),op(map(g->omega(g,aa),G))} minus {aa,op(B)}; #print("A",A):
  B := B union {aa}; #print("B",B):
od;
B;
end:

Bahn_eq := proc(G::list, m, omega, eqproc)
local B, A,aa, g,gaa,i,b;
B := {};
A := {m};
while A <> {} do
  aa := A[1]; #print(op(map(g->omega(g,aa),G)),"next"):
  for g in G do
    gaa := omega(g, aa):
    i := 1;
    while i <= nops(A) and not eqproc(gaa, A[i]) do
      i := i+1;
    od;
    if i > nops(A) then
      A := {op(A), gaa};
    fi;
  od;
  A := A minus {aa};
  for b in B do
    i := 1;
    while i <= nops(A) and not eqproc(b, A[i]) do
      i := i+1;
    od;
    if i <= nops(A) then
      A := subsop(i=NULL, A):
    fi;
  od;
  B := B union {aa}; #print("B",B):
od;
B;
end:


OpPerMa := proc(P::list, M::listlist)
convert(LinearAlgebra[SubMatrix](M, P, P), listlist):
end:

eq_evala := proc(M1, M2) 
evalb({op(map(evala, map(op, convert(M1 - M2, listlist))))} = {0}): 
end proc:


eq_simplify := proc(M1, M2) 
evalb({op(map(simplify, map(op, convert(M1 - M2, listlist))))} = {0}): 
end proc:


h1 := [12, 10, 9, 11, 7, 8, 5, 6, 3, 2, 4, 1];
h2 := [ 1, 6, 3, 5, 4, 2, 11, 10, 9, 8, 7, 12 ]; 
h3 := [ 9, 8, 12, 4, 5, 10, 7, 2, 1, 6, 11, 3 ]; 
h4:= [ 4, 1, 5, 2, 8, 9, 6, 3, 7, 12, 10, 11 ];
h5:=[ 2, 7, 3, 8, 4, 1, 11, 12, 9, 5, 6, 10 ];


s := NewSurface():
DefineEmbedding(s, ddd[1], "faces" = ddd[2], "vertices" = ddd[3]):
coor := CoordinateMatrix(s, 1, "listlist" = true, "radical" = true):
bary := Barycenter(s, 1):
M := Matrix(map(i -> <i - bary>, coor)):
gram := evala((Transpose(M)) . M):
gram_list:=convert(gram,listlist);

AutConditions := proc(g)
local diff, conds, i, j;
#Entrywise difference of Gram matrices
diff := simplify(OpPerMa(g, gram_list) - gram_list);
#Collect nonzero entry conditions
conds := {};
for i to nops(diff) do
    for j to nops(diff[i]) do
        if diff[i][j] <> 0 then
            conds := conds union { diff[i][j]};
        fi;
    od;
od;
return conds;
end:

with(RegularChains):
with(ChainTools): 	
with(SemiAlgebraicSetTools):
R:=PolynomialRing([a,b,c,m,n]);

F:=[op(map(numer,AutConditions(h5))),a>0,b>0,c>0,m<>0,n<>0]; #to find values of a,b,c,m,n for permutation h5
cad:=CylindricalAlgebraicDecompose(F,R,output='rootof');

