!BOI
! !TITLE: Group theoretical determination of the irreducible representations $\Gamma^{\mathcal{G}_0({\bf k})}_{\alpha}$ of the electron eigenstates
! !AUTHORS: Clas Persson
! !AFFILIATION: Royal Institute of Technology, SE--100 44 Stockholm, Sweden
! !INTRODUCTION: Introduction
!   The irreducible representation (IR) of an electron state in a solid denotes the symmetry 
! of the eigenfunction. This symmetry is a fundamental physical quantity of the 
! eigenfunction, which is used -- directly or indirectly -- in various analyses of physical 
! properties. For instance, the transition probability for certain scattering processes can 
! be directly analyzed by knowing the IRs of the eigenstates. 
! Also, understanding the eigenstate degeneracies, and the splitting of degeneneracies due to 
! symmetry-breaking perturbation are direcly connected to the IRs.  
! Determining the IR is useful, not only for analyzing physical properties of the eigenstates, 
! but also when developing and implementing new perturbation Hamiltonians into the computer codes. 
! In these code developments, the determination of the IRs serves as a tool for insuring that 
! the symmetries of the eigenfunctions are correct after the new implementation; this is a 
! fairly sensitive test. 
! Moreover, since it is expected that two electronic energy bands with same IR do not 
! "intersected", the IRs can be utilized to predict and display smoother energy bands with less 
! number of {\bf k}-points along the symmetry lines. This saves computational time, as well as 
! saves disc space for large systems.  
!
! \vspace{24pt}
! The code implementation as well as the notation of the IRs is according to:  
!
! C. Persson, Comput. Phys. Commun. {\bf XX}, XX (2008),  
!
! wherein properties of the IRs are illustrated and references are given. 
! In order to present well-defined IRs, it is important to know that:
! \vspace{-8pt}
! \begin{itemize}
!   \item Labeling of IRs depends on choise of crystal origin and orientation. 
!   \item surface
!   \item 
!   \item 
! \end{itemize}
!  
! 
!    \vspace{24pt}
!    Clas Persson
!
!    \vspace{12pt}
!    Stockholm, September 2007
!    \newpage
!
! \section{Basic Definitions}
! Below, a very brief introduction to the theory of IRs is presented, and 
! the most important definitions are given. A somewhat more complete description 
! of the theories are found in C. Persson, Comput. Phys. Commun. {\bf XX}, XX (2008), 
! originally following the discussion by 
! J. F. Cornwell, "Group Theory in Physics", Vol. 1 (Academic Press, London, 1984). 
!
! \subsection{Crystallographic Space Groups}
! The set of all coordinate transformations ${\bf P}_g$ ($g$ = 1, 2, ... $\infty$) that maps 
! the atomic positions in an infinite crystal into itself forms an infinite 
! {\it crystallographic space group} $\mathcal{G}_{\infty}$. 
! An element ${\bf P}_g$ = $\{{\bf R}_g|{\bf t}_g\}$ of the space group consists of a
! rotation part ${\bf R}_g$ and a translation part ${\bf t}_g$. 
! By definition, $[{\bf P}_g, \hat{H}]$ = 0.
! The inverse transformation ${\bf P}_g^{-1}$ equals $\{{\bf R}_g^{-1}|-{\bf R}_g^{-1}{\bf t}_g\}$.
! Operating $\{{\bf R}_g|{\bf t}_g\}$ $\in$ $\mathcal{G}_{\infty}$ 
! on an eigenfunction $\psi_i({\bf r})$ to the Hamiltonian of the crystal implies: 
! 
! \begin{equation}\label{e:P_operate}
! {\bf P}_g\ \psi_i({\bf r}) = \{{\bf R}_g|{\bf t}_g\}\psi_i({\bf r}) = 
! \psi_i({\{{\bf R}_g|{\bf t}_g\}^{-1}\bf r}) = 
! \psi_i({\bf R}_g^{-1}{\bf r}-{\bf R}_g^{-1}{\bf t}_g).
! \end{equation}
! 
! The set of all pure rotations $\{{\bf R}_g|{\bf 0}\}$ forms a finite {\it crystallographic point group} 
! $\mathcal{G}_0$ of the space group $\mathcal{G}_{\infty}$. 
! To every coordinate transformation $\{{\bf R}_g|{\bf 0}\}$  $\in$ $\mathcal{G}_0$ there 
! exists a unique (up to a lattice vector) non-primitive vector ${\bf \tau}_g$ such that 
! $\{{\bf R}_g|{\bf t}_g\}$ = $\{{\bf R}_g|{\bf \tau}_g+{\bf t}\}$ $\in$ $\mathcal{G}_{\infty}$ 
! for every crystallographic lattice vector {\bf t} = $n_1{\bf a}_1$ + $n_2{\bf a}_2$ + $n_3{\bf a}_3$, 
! where $n_\iota$ are integer numbers and ${\bf a}_\iota$ are the 
! primitive translation vectors of the crystallographic lattice.
! For definiteness, we choose ${\bf \tau}_g$ to be the smallest possible 
! non-primitive vector: ${\bf \tau}_g$ = 
! $q_1{\bf a}_1$ + $q_2{\bf a}_2$ + $q_3{\bf a}_3$ with 0 $\le$  $q_\iota$ $<$ 1.
! 
!
! The set of all pure primitive translations $\{{\bf 1}|{\bf t}\}$ forms an infinite group 
! $\mathcal{T}_{\infty}$. 
! However, the Born--von Karman periodic boundary condition of a lattice 
! reduces the infinite translation group $\mathcal{T}_{\infty}$ to a finite group $\mathcal{T}$ of
! pure primitive translations $\{{\bf 1}|{\bf t}_n\}$ with the order of $N$. 
! The corresponding finite space group $\mathcal{G}$ contains the 
! transformations $\{{\bf R}_g|{\bf \tau}_g+{\bf t}_n\}$ and has the same point group
! $\mathcal{G}_0$ as $\mathcal{G}_{\infty}$. 
! Moreover, in many cases (see Section \ref{s:ALLOW}) the set of coordinate transformations 
! $\{{\bf R}_g|{\bf \tau}_g+{\bf t}_n\}$ $\in$ $\mathcal{G}$ can be further reduced to the set of 
! transformations $\{{\bf R}_g|{\bf \tau}_g\}$ $\in$ $\mathcal{G}$, with the same number of 
! symmetry transformations as number of different pure rotations.
!
! A {\it class} is generated by ${\bf P}_{g'}$${\bf P}_{g}$${\bf P}_{g'}^{-1}$, 
! where $g'$ goes through all the group elements of $\mathcal{G}$. 
!
! \subsection{Irreducible Representations}\label{s:IR}
! Consider $d$ linearly independent functions $\psi_1({\bf r})$, 
! $\psi_2({\bf r})$, ... $\psi_d({\bf r})$, spanning the particle space $W$. 
! If all elements  ${\bf P}_g$ $\in$ $\mathcal{G}$ operating on $\psi_i({\bf r})$ 
! form functions which all lie in $W$, i.e.,
! 
! \begin{equation}\label{e:P_operateG}
! {\bf P}_g\ \psi_i({\bf r}) = \sum_{i'} \Gamma^{\mathcal{G}}({\bf P}_g)_{i'i}\ \psi_{i'}({\bf r})  
! \ \ \ \ \ \ \ \ \ \ \ \ i,i' = 1, 2..., d,
! \end{equation}
! 
! \noindent then $\psi_i({\bf r})$ form a basis for space $W$. 
! Each element ${\bf P}_g$ can be assigned to a non-singular 
! $d$$\times$$d$ matrix $\Gamma^{\mathcal{G}}({\bf P}_g)$ which represents the group element 
! ${\bf P}_g$ in $W$.
! The set of matrices (i.e., the representation matrices of all ${\bf P}_g$) forms the 
! $d$-dimensional {\it representation} $\Gamma^{\mathcal{G}}$ of $\mathcal{G}$.
! The subspace $W$ is {\it reducible} if there exists a subspace $W_{\alpha}$ of $W$, invariant 
! under the group $\mathcal{G}$ and spanned by the $d_{\alpha}$ linearly independent basis functions 
! $\psi_{i_{\alpha}+1}({\bf r})$, $\psi_{i_{\alpha}+2}({\bf r})$, ... $\psi_{i_{\alpha}+d_{\alpha}}({\bf r})$. 
! If the subspace is completely reducible one can always find a block diagonalized form of 
! $\Gamma^{\mathcal{G}}({\bf P}_g)$ for every ${\bf P}_g$, containing representation matrices
! $\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)$ that cannot be further block diagonalized:
! 
! \begin{equation}\label{e:G_blockAB2}
! {\bf P}_g\ \psi_{i_{\alpha}+i}({\bf r}) = \sum_{i'} \Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)_{i'i}
!     \ \psi_{i_{\alpha}+i'}({\bf r}) \ \ \ \ \ \ \ \ \ \ \ \ \left\{
!          \begin{array}[]{lr}
!                 i,i'  & = 1,2,...,d_{\alpha} \\
!                \alpha & = 1,2,...,\alpha_n \\
!          \end{array} 
!          \right.
! \end{equation}
! \begin{equation}\label{e:G_blockAB}
! \Gamma^{\mathcal{G}}({\bf P}_g) = \left(
!   \begin{array}[]{cccc}
! 	  \Gamma^{\mathcal{G}}_{1}({\bf P}_g) & {\bf 0} & \vdots & {\bf 0} \\
!   	{\bf 0} & \Gamma^{\mathcal{G}}_{2}({\bf P}_g) & \vdots & {\bf 0} \\
!   	\cdots & \cdots & \vdots & \cdots \\
! 	  {\bf 0} & {\bf 0} & \vdots & \: \ \Gamma^{\mathcal{G}}_{\alpha_n}({\bf P}_g) 
!   \end{array} 
! \right) \\  
! \end{equation}
! 
! 
! \noindent The $d_{\alpha}$$\times$$d_{\alpha}$ matrix $\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)$ 
! representing the group element ${\bf P}_g$ in $W_{\alpha}$ is not unique; it depends on the choice 
! of basis functions. However, the trace of a representation matrix (i.e., the {\it character}) of 
! $\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)$ is independent on the choice of basis functions:
! 
! \begin{equation}\label{e:G_trace}
!   \chi^{\mathcal{G}}_{\alpha}({\bf P}_g) = tr\left[\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)\right] = 
!  \sum_{i} \Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)_{ii}. 
! \end{equation}
! 
! For a given $\alpha$, the set of $d_{\alpha}$$\times$$d_{\alpha}$ representation matrices 
! $\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)$ for all ${\bf P}_g$, forms the
! $d_{\alpha}$-dimensional IR $\Gamma^{\mathcal{G}}_{\alpha}$ of the space group $\mathcal{G}$. 
! The characters $\chi^{\mathcal{G}}_{\alpha}({\bf P}_g)$ are a signature for the symmetry of the 
! basis functions 
! $\psi_{i_{\alpha}+1}({\bf r})$, $\psi_{i_{\alpha}+2}({\bf r})$, ... $\psi_{i_{\alpha}+d_{\alpha}}({\bf r})$, 
! and how these functions transform under the operation ${\bf P}_g$. 
! All representation matrices within a class have equal character, and therefore ${\bf P}_g$ are
! grouped into their classes.
!
! The reducible representation $\Gamma^{\mathcal{G}}$ of $\mathcal{G}$ can be expressed in terms of 
! the IRs by $\Gamma^{\mathcal{G}}$ = $\Gamma^{\mathcal{G}}_1$$\:\oplus\:$$\Gamma^{\mathcal{G}}_2$
! $\:\oplus$, ...,$\oplus\:$$\Gamma^{\mathcal{G}}_{\alpha_n}$, where the {\it direct sum} operation 
! $\oplus$ is defined by the block-diagonal form in Eq.~(\ref{e:G_blockAB}).
! If $\Gamma^{\mathcal{G}}_{\alpha}$ and $\Gamma^{\mathcal{G}}_{\alpha'}$ are two 
! IRs of the group $\mathcal{G}$, having the $d_{\alpha}$$\times$$d_{\alpha}$ and
! $d_{\alpha'}$$\times$$d_{\alpha'}$ matrices $\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)$
! and $\Gamma^{\mathcal{G}}_{\alpha'}({\bf P}_g)$, respectively, then the {\it direct} 
! {\it product representation} $\Gamma^{\mathcal{G}}_{\alpha}$$\otimes$$\Gamma^{\mathcal{G}}_{\alpha'}$ is
! defined via the $d_{\alpha}$$d_{\alpha'}$$\times$$d_{\alpha}$$d_{\alpha'}$ representation matrices
! $\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)$$\otimes$$\Gamma^{\mathcal{G}}_{\alpha'}({\bf P}_g)$ with the 
! Kronecker (or tensor) matrix elements
!  $\left[\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)\otimes\Gamma^{\mathcal{G}}_{\alpha'}({\bf P}_g)\right]_
!  {(i'j')(ij)}$ = $\Gamma^{\mathcal{G}}_{\alpha}({\bf P}_g)_{i'i}$ $\cdot$ 
!                 $\Gamma^{\mathcal{G}}_{\alpha'}({\bf P}_g)_{j'j}$.
!
!
! \subsection{Single and Double Groups}\label{s:DG}  
! If the Hamiltonian involves interactions through the spin component of the electrons, 
! the basis functions have to be represented by vector functions (or spinors). 
! Here it is sufficient to use two-component spinors used in the present code:
! 
! \begin{equation}\label{e:DG_spinor}
! \bar{\psi}_{j_{\alpha}{\bf k}}({\bf r}) =  \left(
!        \begin{array}[]{l}
!          \psi^{\uparrow}_{j_{\alpha}{\bf k}}({\bf r}) \\ 
!          \psi^{\downarrow}_{j_{\alpha}{\bf k}}({\bf r}) \\
!        \end{array} 
! \right) \ \ \ \ j_{\alpha} = l_{\alpha}+1, ..., l_{\alpha}+d_{\alpha} 
! \end{equation}
! 
! When a coordinate transformation ${\bf P}_g$ = $\{{\bf R}_g|{\bf \tau}_g + {\bf t}_n\}$ is 
! operating on a spinor in the real space, also a transformation ${\bf U}({\bf R}_g^p)$ is 
! operating on spin components in the spin space
! 
! \begin{eqnarray}\label{e:DG_operate}
!  \pm{\bf U}({\bf R}_g^p){\bf P}_g\ \bar{\psi}_{j_{\alpha}{\bf k}}({\bf r}) =
!  \bar{\psi}'_{j_{\alpha}{\bf k}}({\bf r}') = 
!  \sum_{j'} \Gamma^{\mathcal{G}^D}_{\alpha}({\bf P}_g)_{jj'}
! \bar{\psi}_{j'_{\alpha}{\bf k}}({\bf r}) \  \ \ \left\{
!        \begin{array}[]{ll}
!          j_{\alpha} & = l_{\alpha}+1, ..., l_{\alpha}+d_{\alpha} \\ 
!          \alpha     & = 1,... \alpha_n.  
!        \end{array} 
! \right.
! \end{eqnarray}
! 
! where ${\bf R}_g^p$ = $|{\bf R}_g| \cdot {\bf R}_g$
! is the {\it proper rotation} of ${\bf R}_g$.
! The sign $\pm$ is entirely a matter of convention, but it 
! implies twice as many operators as for case of the spin-independent 
! eigenfunctions, and the space group $\mathcal{G}^D$ corresponding to spinors is therefore 
! called a {\it double group} compared to the {\it single group} $\mathcal{G}$ for  
! spin-independent eigenfunctions.  
! The 2$\times$2 spin-operation matrices ${\bf U}({\bf R}_g^p)$ are related to ${\bf R}_g^p$
! through the Euler's angles $\varphi$, $\theta$, and $\xi$ as:
!
! \begin{eqnarray}\label{e:DG_euler}
! {\bf U}({\bf R}_g^p) = && \left(
!       \begin{array}[]{rr}
!         cos(\theta/2)e^{+i(\xi+\varphi)/2} & sin(\theta/2)e^{+i(\xi-\varphi)/2}\\ 
!        -sin(\theta/2)e^{-i(\xi-\varphi)/2} & cos(\theta/2)e^{-i(\xi+\varphi)/2} \\
!       \end{array} 
! \right) \\ \nonumber && \\ \nonumber
! && 0\leq\theta\leq\pi, \ 0\leq\xi\leq4\pi, \ \textnormal{ and } 0\leq\varphi\leq\pi.
! \end{eqnarray}
!
! \noindent where we follow the following convention: first a rotation trough the angle $\varphi$ 
! about the $z$-axis, then a rotation through the angle $\theta$ about the new $y$-axis, 
! and last a rotation through the angle $\xi$ about the new $z$-axis.
!
! \subsection{The Group of Allowed Wave Vector {\bf k}}\label{s:ALLOW} 
! For all coordinate transformations ${\bf P}_g$ = $\{{\bf R}_g|{\bf t}_g\}$ $\in$ $\mathcal{G}$,
! the pure rotation operations $\{{\bf R}_g|{\bf 0}\}$ $\in$ $\mathcal{G}_0$
! possess the property that for any reciprocal lattice vector ${\bf K}$, the vector
! ${\bf R}_g{\bf K}$ is a lattice vector ${\bf K}'$ of the same reciprocal lattice.
! The {\it group of allowed wave vector {\bf k}} has the 
! property ${\bf R}_g{\bf k}$ = ${\bf k}$ + ${\bf K}$. 
! An IR of the space group $\mathcal{G}({\bf k})$ can in most cases be expressed in terms of 
! an IR of the corresponding point group $\mathcal{G}_0({\bf k})$. If
! $ e^{-i{\bf k} ( {\bf R}_g{\bf \tau}_{g'} - {\bf \tau}_{g'})}$ = 1
!  is true for every pair $\{ \{{\bf R}_g|{\bf \tau}_g\};\{{\bf R}_{g'}|{\bf \tau}_{g'}\} \}$ 
! $\in$ $\mathcal{G}({\bf k})$, then 
! 
! \begin{equation}\label{e:Allowed2}
! \Gamma^{\mathcal{G}({\bf k})}_{\alpha}(\{{\bf R}_g|{\bf \tau}_g+{\bf t}_n\})   = 
! e^{-i{\bf k} ( {\bf \tau}_g + {\bf t}_n )} \cdot
! \Gamma^{\mathcal{G}_0({\bf k})}_{\alpha}(\{{\bf R}_g|{\bf 0}\}), 
! \end{equation}
! 
! $\Gamma^{\mathcal{G}_0({\bf k})}_{\alpha}$ is therefore called the {\it relevant} IR 
! of $\mathcal{G}({\bf k})$.
! This is a very useful relation between the space-group IR and the point-group IR, commonly 
! used when presenting the IRs in electronic band structures.
! The relation of Eq.~\ref{e:Allowed2} holds for all symmorphic space groups, and 
! also for every internal point of the Brillouin zone. However, it can (but it does not have to)
! be invalid on the surface of the  for non-symmorphic space groups. 
! In those cases there is an extended method based on factor groups, which to date is {\it not} 
! implemented in the code (since it is not trivial to define labels for those IRs).
! Presenting the relevant IRs, without the phase $exp(-i{\bf k}({\bf \tau}_g+{\bf t}_n ))$,
! causes sometimes an ambiguous labeling of the IRs, depending on crystal origin and orientation.
! The present code therefore ouputs also the phase for each class.
! 
! \subsection{Time Inversion}\label{s:TR} 
! Physical systems should be {\it invariant under time inversion}. This symmetry is in general 
! not comprised in the set of coordinate transformations, and must thereby be considered
! additionally by the time-reversal operator $\hat{K}$.
! For spin-independent systems $\hat{K}$ = $\hat{K}_0$, where $\hat{K}_0$ is defined by 
! $\hat{K}_0\psi_{i}({\bf r})$ = $\psi^*_{i}({\bf r})$
! For spin-dependent systems $\hat{K}$ can be represented by $\sigma_2\hat{K}_0$
! where $\sigma$ is the vector of the three Pauli spin matrices:
!
! \begin{equation}\label{e:TR_spinor}
! \hat{K}\bar{\psi}_{j_{\alpha}{\bf k}}({\bf r}) = \sigma_2\hat{K}_0 \left(
!        \begin{array}[]{r}
!         \psi^{\uparrow}_{j_{\alpha}{\bf k}}({\bf r}) \\
!          \psi^{\downarrow}_{j_{\alpha}{\bf k}}({\bf r})
!        \end{array} 
!        \right) = \left(
!        \begin{array}[]{r}
!          -i\psi^{\downarrow *}_{j_{\alpha}{\bf k}}({\bf r}) \\
!           i\psi^{\uparrow *}_{j_{\alpha}{\bf k}}({\bf r})
!        \end{array} 
!        \right) = \bar{\psi}'_{j_{\alpha}{\bf k}}({\bf r}).
! \end{equation}
! 
! Either 
! $\bar{\psi}'_{j_{\alpha}{\bf k}}({\bf r})$ are linear combinations of 
! $\bar{\psi}_{j_{\alpha}{\bf k}}({\bf r})$ and thereby span the same subspace $W^D_{\alpha}({\bf k})$, or 
! $\bar{\psi}'_{j_{\alpha}{\bf k}}({\bf r})$ are linearly independent of 
! $\bar{\psi}_{j_{\alpha}{\bf k}}({\bf r})$ and thereby span another subspace $W^D_{\alpha'}({\bf k})$.
! In the latter case,the two different subspaces $W^D_{\alpha}({\bf k})$ and $W^D_{\alpha'}({\bf k})$ 
! are degenerate due to the time-reversal symmetry.
!
!
!EOI

