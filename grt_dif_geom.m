(* ::Package:: *)

(* First, we introduce the coordinates. *)
Dim = 4; (* Total number of dimensions, only intended to be equal 4, but who knows... *)
crd = {t,r,\[Theta],\[CurlyPhi]}; (* List of the coordinate labels. The number of coordinates must match the dimension, obviously. *)


(* Here we introduce the metric tensor. *)

gCov = Table[0,{Dim},{Dim}]; (* Covariant metric tensor. *)
gCon = Table[0,{Dim},{Dim}]; (* Contravariant metric tensor. *)

(* The only required input is the metric components. *)
(* For simplicity the indices run from 1 to 4, and 1 corresponds to temporal dimension *)
gCov[[1,1]] = 1-a/crd[[2]];
gCov[[2,2]] = -(1/(1-a/crd[[2]]));
gCov[[3,3]]= - crd[[2]]^2;
gCov[[4,4]]= - ( crd[[2]] * Sin[ crd[[3]] ] ) ^ 2;


(* Contravariant metric is inverse of the covariant. *)
gCon = Simplify[Inverse[gCov]];


(* Our very own Kronecker delta, just in case. *)

delta = Table[0,{Dim},{Dim}];
Do[delta[[iind,iind]] = 1, {iind, Dim}];
delta


(* Here we declare various geomentic objects.*)

Christoffel = Table[0,{Dim},{Dim},{Dim}];   (* Christoffel symbols *)
Riemann = Table[0,{Dim},{Dim},{Dim},{Dim}]; (* Curvature tensor *)
Ricci = Table[0,{Dim},{Dim}]; (* Ricci tensor *)
Einstein = Table[0,{Dim},{Dim}]; (* Einstein tensor *)
ScalarCurv = 0; (*Scalar curvature*)
Bianchi = Table[0,{Dim}]; (*Bianchi identities *)

(* 
With regeards to the placement of the indices and sign conventions, we use the following definitions: 

Christoffel = \Gamma ^i_{jk}
Riemann = Riemann^i_{jkl}
Ricci = Ricci_ {ij} = Riemann^m_{imj}
ScalarCurv = gCov^{ij} Ricci_ {ij}
Einstein = gCon^{im} Ricci_ {jm} - 1/2 delta^i_j  ScalarCurv
Bianchi = Einstein^i_{j;i}
*)


Christoffel = Simplify[
	Table[
		1/2 * Sum[
				gCon[[iind,nind]] * (
					D[ gCov[[jind,nind]] , crd[[kind]] ]
					+ D[ gCov[[kind,nind]] , crd[[jind]] ]
					- D[ gCov[[jind,kind]] , crd[[nind]] ] 
				)
		, {nind, 1, Dim} ]
	, {iind, 1, Dim}, {jind, 1, Dim}, {kind, 1, Dim} ]
]


Riemann = Simplify[
		Table[
			D[ Christoffel[[iind,jind,lind]] , crd[[kind]] ]
			- D[ Christoffel[[iind,jind,kind]] , crd[[lind]] ]
			+ Sum[ Christoffel[[iind,nind,kind]] * Christoffel[[nind,jind,lind]] , {nind, 1, Dim} ]
			- Sum[ Christoffel[[iind,nind,lind]] * Christoffel[[nind,jind,kind]] , {nind, 1, Dim} ]
		, {iind, 1, Dim}, {jind, 1, Dim}, {kind, 1, Dim}, {lind, 1, Dim} ]
]



Ricci = Simplify[
	Table[
		Sum[ Riemann[[nind,jind,nind,lind]] , {nind, 1, Dim} ] 
	, {jind , 1, Dim}, {lind, 1, Dim} ] 
]


ScalarCurv = Simplify[
	Sum[ Sum[ 
		gCon[[mind,nind]] * Ricci[[mind,nind]] 
	, {nind, 1, Dim} ], {mind, 1, Dim} ]
]


Einstein = Simplify[
	Table[
		Sum[ gCon[[iind,nind]] * Ricci[[jind,nind]] , {nind, 1, Dim} ]
		- 1/2 * delta[[iind,jind]] * ScalarCurv
		
	, {iind , 1, Dim}, {jind, 1, Dim} ]
]


Bianchi = Simplify[
	Table[
		Sum[ D[ Einstein[[lind,iind]] , crd[[lind]] ] , {lind, 1, Dim} ]
		+ Sum[ Christoffel[[lind,mind,lind]] * Einstein[[mind,iind]] , {lind, 1, Dim}, {mind, 1, Dim} ]
		- Sum[ Christoffel[[mind,iind,lind]] * Einstein[[lind,mind]] , {lind, 1, Dim}, {mind, 1, Dim} ]
	, {iind, 1, Dim} ]
]



Clear[crd]; Clear[gCov]; Clear[gCon]; Clear[Christoffel]; Clear[Riemann]; Clear[Ricci]; Clear[ScalarCurv]; 
Clear[Einstein]; Clear[Bianchi]; Clear[delta]; Clear[Dim];

