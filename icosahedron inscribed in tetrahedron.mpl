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



