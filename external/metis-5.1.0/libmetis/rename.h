/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 * $Id: rename.h 13933 2013-03-29 22:20:46Z karypis $
 *
 */


#ifndef _LIBMETIS_RENAME_H_
#define _LIBMETIS_RENAME_H_


/* balance.c */
#define Balance2Way			libmetis_MUMPS_Balance2Way
#define Bnd2WayBalance			libmetis_MUMPS_Bnd2WayBalance
#define General2WayBalance		libmetis_MUMPS_General2WayBalance
#define McGeneral2WayBalance            libmetis_MUMPS_McGeneral2WayBalance

/* bucketsort.c */
#define BucketSortKeysInc		libmetis_MUMPS_BucketSortKeysInc

/* checkgraph.c */
#define CheckGraph                      libmetis_MUMPS_CheckGraph
#define CheckInputGraphWeights          libmetis_MUMPS_CheckInputGraphWeights
#define FixGraph                        libmetis_MUMPS_FixGraph

/* coarsen.c */
#define CoarsenGraph			libmetis_MUMPS_CoarsenGraph
#define Match_RM                        libmetis_MUMPS_Match_RM
#define Match_SHEM                      libmetis_MUMPS_Match_SHEM
#define Match_2Hop                      libmetis_MUMPS_Match_2Hop
#define Match_2HopAny                   libmetis_MUMPS_Match_2HopAny
#define Match_2HopAll                   libmetis_MUMPS_Match_2HopAll
#define PrintCGraphStats                libmetis_MUMPS_PrintCGraphStats
#define CreateCoarseGraph		libmetis_MUMPS_CreateCoarseGraph
#define CreateCoarseGraphNoMask		libmetis_MUMPS_CreateCoarseGraphNoMask
#define CreateCoarseGraphPerm		libmetis_MUMPS_CreateCoarseGraphPerm
#define SetupCoarseGraph		libmetis_MUMPS_SetupCoarseGraph
#define ReAdjustMemory			libmetis_MUMPS_ReAdjustMemory

/* compress.c */
#define CompressGraph			libmetis_MUMPS_CompressGraph
#define PruneGraph			libmetis_MUMPS_PruneGraph

/* contig.c */
#define FindPartitionInducedComponents  libmetis_MUMPS_FindPartitionInducedComponents   
#define IsConnected                     libmetis_MUMPS_IsConnected
#define IsConnectedSubdomain            libmetis_MUMPS_IsConnectedSubdomain
#define FindSepInducedComponents        libmetis_MUMPS_FindSepInducedComponents
#define EliminateComponents             libmetis_MUMPS_EliminateComponents
#define MoveGroupContigForCut           libmetis_MUMPS_MoveGroupContigForCut
#define MoveGroupContigForVol           libmetis_MUMPS_MoveGroupContigForVol

/* debug.c */
#define ComputeCut			libmetis_MUMPS_ComputeCut
#define ComputeVolume			libmetis_MUMPS_ComputeVolume
#define ComputeMaxCut			libmetis_MUMPS_ComputeMaxCut
#define CheckBnd			libmetis_MUMPS_CheckBnd
#define CheckBnd2			libmetis_MUMPS_CheckBnd2
#define CheckNodeBnd			libmetis_MUMPS_CheckNodeBnd
#define CheckRInfo			libmetis_MUMPS_CheckRInfo
#define CheckNodePartitionParams	libmetis_MUMPS_CheckNodePartitionParams
#define IsSeparable			libmetis_MUMPS_IsSeparable
#define CheckKWayVolPartitionParams     libmetis_MUMPS_CheckKWayVolPartitionParams

/* fm.c */
#define FM_2WayRefine                   libmetis_MUMPS_FM_2WayRefine
#define FM_2WayCutRefine                libmetis_MUMPS_FM_2WayCutRefine
#define FM_Mc2WayCutRefine              libmetis_MUMPS_FM_Mc2WayCutRefine
#define SelectQueue                     libmetis_MUMPS_SelectQueue
#define Print2WayRefineStats            libmetis_MUMPS_Print2WayRefineStats

/* fortran.c */
#define Change2CNumbering		libmetis_MUMPS_Change2CNumbering
#define Change2FNumbering		libmetis_MUMPS_Change2FNumbering
#define Change2FNumbering2		libmetis_MUMPS_Change2FNumbering2
#define Change2FNumberingOrder		libmetis_MUMPS_Change2FNumberingOrder
#define ChangeMesh2CNumbering		libmetis_MUMPS_ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		libmetis_MUMPS_ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		libmetis_MUMPS_ChangeMesh2FNumbering2

/* graph.c */
#define SetupGraph			libmetis_MUMPS_SetupGraph
#define SetupGraph_adjrsum              libmetis_MUMPS_SetupGraph_adjrsum
#define SetupGraph_tvwgt                libmetis_MUMPS_SetupGraph_tvwgt
#define SetupGraph_label                libmetis_MUMPS_SetupGraph_label
#define SetupSplitGraph                 libmetis_MUMPS_SetupSplitGraph
#define CreateGraph                     libmetis_MUMPS_CreateGraph
#define InitGraph                       libmetis_MUMPS_InitGraph
#define FreeRData                       libmetis_MUMPS_FreeRData
#define FreeGraph                       libmetis_MUMPS_FreeGraph

/* initpart.c */
#define Init2WayPartition		libmetis_MUMPS_Init2WayPartition
#define InitSeparator			libmetis_MUMPS_InitSeparator
#define RandomBisection			libmetis_MUMPS_RandomBisection
#define GrowBisection			libmetis_MUMPS_GrowBisection
#define McRandomBisection               libmetis_MUMPS_McRandomBisection
#define McGrowBisection                 libmetis_MUMPS_McGrowBisection
#define GrowBisectionNode		libmetis_MUMPS_GrowBisectionNode

/* kmetis.c */
#define MlevelKWayPartitioning		libmetis_MUMPS_MlevelKWayPartitioning
#define InitKWayPartitioning            libmetis_MUMPS_InitKWayPartitioning

/* kwayfm.c */
#define Greedy_KWayOptimize		libmetis_MUMPS_Greedy_KWayOptimize
#define Greedy_KWayCutOptimize		libmetis_MUMPS_Greedy_KWayCutOptimize
#define Greedy_KWayVolOptimize          libmetis_MUMPS_Greedy_KWayVolOptimize
#define Greedy_McKWayCutOptimize        libmetis_MUMPS_Greedy_McKWayCutOptimize
#define Greedy_McKWayVolOptimize        libmetis_MUMPS_Greedy_McKWayVolOptimize
#define IsArticulationNode              libmetis_MUMPS_IsArticulationNode
#define KWayVolUpdate                   libmetis_MUMPS_KWayVolUpdate

/* kwayrefine.c */
#define RefineKWay			libmetis_MUMPS_RefineKWay
#define AllocateKWayPartitionMemory	libmetis_MUMPS_AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	libmetis_MUMPS_ComputeKWayPartitionParams
#define ProjectKWayPartition		libmetis_MUMPS_ProjectKWayPartition
#define ComputeKWayBoundary		libmetis_MUMPS_ComputeKWayBoundary
#define ComputeKWayVolGains             libmetis_MUMPS_ComputeKWayVolGains
#define IsBalanced			libmetis_MUMPS_IsBalanced

/* mcutil */
#define rvecle                          libmetis_MUMPS_rvecle
#define rvecge                          libmetis_MUMPS_rvecge
#define rvecsumle                       libmetis_MUMPS_rvecsumle
#define rvecmaxdiff                     libmetis_MUMPS_rvecmaxdiff
#define ivecle                          libmetis_MUMPS_ivecle
#define ivecge                          libmetis_MUMPS_ivecge
#define ivecaxpylez                     libmetis_MUMPS_ivecaxpylez
#define ivecaxpygez                     libmetis_MUMPS_ivecaxpygez
#define BetterVBalance                  libmetis_MUMPS_BetterVBalance
#define BetterBalance2Way               libmetis_MUMPS_BetterBalance2Way
#define BetterBalanceKWay               libmetis_MUMPS_BetterBalanceKWay
#define ComputeLoadImbalance            libmetis_MUMPS_ComputeLoadImbalance
#define ComputeLoadImbalanceDiff        libmetis_MUMPS_ComputeLoadImbalanceDiff
#define ComputeLoadImbalanceDiffVec     libmetis_MUMPS_ComputeLoadImbalanceDiffVec
#define ComputeLoadImbalanceVec         libmetis_MUMPS_ComputeLoadImbalanceVec

/* mesh.c */
#define CreateGraphDual                 libmetis_MUMPS_CreateGraphDual
#define FindCommonElements              libmetis_MUMPS_FindCommonElements
#define CreateGraphNodal                libmetis_MUMPS_CreateGraphNodal
#define FindCommonNodes                 libmetis_MUMPS_FindCommonNodes
#define CreateMesh                      libmetis_MUMPS_CreateMesh
#define InitMesh                        libmetis_MUMPS_InitMesh
#define FreeMesh                        libmetis_MUMPS_FreeMesh

/* meshpart.c */
#define InduceRowPartFromColumnPart     libmetis_MUMPS_InduceRowPartFromColumnPart

/* minconn.c */
#define ComputeSubDomainGraph           libmetis_MUMPS_ComputeSubDomainGraph
#define UpdateEdgeSubDomainGraph        libmetis_MUMPS_UpdateEdgeSubDomainGraph
#define PrintSubDomainGraph             libmetis_MUMPS_PrintSubDomainGraph
#define EliminateSubDomainEdges         libmetis_MUMPS_EliminateSubDomainEdges
#define MoveGroupMinConnForCut          libmetis_MUMPS_MoveGroupMinConnForCut
#define MoveGroupMinConnForVol          libmetis_MUMPS_MoveGroupMinConnForVol

/* mincover.c */
#define MinCover			libmetis_MUMPS_MinCover
#define MinCover_Augment		libmetis_MUMPS_MinCover_Augment
#define MinCover_Decompose		libmetis_MUMPS_MinCover_Decompose
#define MinCover_ColDFS			libmetis_MUMPS_MinCover_ColDFS
#define MinCover_RowDFS			libmetis_MUMPS_MinCover_RowDFS

/* mmd.c */
#define genmmd				libmetis_MUMPS_genmmd
#define mmdelm				libmetis_MUMPS_mmdelm
#define mmdint				libmetis_MUMPS_mmdint
#define mmdnum				libmetis_MUMPS_mmdnum
#define mmdupd				libmetis_MUMPS_mmdupd


/* ometis.c */
#define MlevelNestedDissection		libmetis_MUMPS_MlevelNestedDissection
#define MlevelNestedDissectionCC	libmetis_MUMPS_MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	libmetis_MUMPS_MlevelNodeBisectionMultiple
#define MlevelNodeBisectionL2		libmetis_MUMPS_MlevelNodeBisectionL2
#define MlevelNodeBisectionL1		libmetis_MUMPS_MlevelNodeBisectionL1
#define SplitGraphOrder			libmetis_MUMPS_SplitGraphOrder
#define SplitGraphOrderCC		libmetis_MUMPS_SplitGraphOrderCC
#define MMDOrder			libmetis_MUMPS_MMDOrder

/* options.c */
#define SetupCtrl                       libmetis_MUMPS_SetupCtrl
#define SetupKWayBalMultipliers         libmetis_MUMPS_SetupKWayBalMultipliers
#define Setup2WayBalMultipliers         libmetis_MUMPS_Setup2WayBalMultipliers
#define PrintCtrl                       libmetis_MUMPS_PrintCtrl
#define FreeCtrl                        libmetis_MUMPS_FreeCtrl
#define CheckParams                     libmetis_MUMPS_CheckParams

/* parmetis.c */
#define MlevelNestedDissectionP		libmetis_MUMPS_MlevelNestedDissectionP
#define FM_2WayNodeRefine1SidedP        libmetis_MUMPS_FM_2WayNodeRefine1SidedP
#define FM_2WayNodeRefine2SidedP        libmetis_MUMPS_FM_2WayNodeRefine2SidedP

/* pmetis.c */
#define MlevelRecursiveBisection	libmetis_MUMPS_MlevelRecursiveBisection
#define MultilevelBisect		libmetis_MUMPS_MultilevelBisect
#define SplitGraphPart			libmetis_MUMPS_SplitGraphPart

/* refine.c */
#define Refine2Way			libmetis_MUMPS_Refine2Way
#define Allocate2WayPartitionMemory	libmetis_MUMPS_Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	libmetis_MUMPS_Compute2WayPartitionParams
#define Project2WayPartition		libmetis_MUMPS_Project2WayPartition

/* separator.c */
#define ConstructSeparator		libmetis_MUMPS_ConstructSeparator
#define ConstructMinCoverSeparator	libmetis_MUMPS_ConstructMinCoverSeparator

/* sfm.c */
#define FM_2WayNodeRefine2Sided         libmetis_MUMPS_FM_2WayNodeRefine2Sided 
#define FM_2WayNodeRefine1Sided         libmetis_MUMPS_FM_2WayNodeRefine1Sided
#define FM_2WayNodeBalance              libmetis_MUMPS_FM_2WayNodeBalance

/* srefine.c */
#define Refine2WayNode			libmetis_MUMPS_Refine2WayNode
#define Allocate2WayNodePartitionMemory	libmetis_MUMPS_Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	libmetis_MUMPS_Compute2WayNodePartitionParams
#define Project2WayNodePartition	libmetis_MUMPS_Project2WayNodePartition

/* stat.c */
#define ComputePartitionInfoBipartite   libmetis_MUMPS_ComputePartitionInfoBipartite
#define ComputePartitionBalance		libmetis_MUMPS_ComputePartitionBalance
#define ComputeElementBalance		libmetis_MUMPS_ComputeElementBalance

/* timing.c */
#define InitTimers			libmetis_MUMPS_InitTimers
#define PrintTimers			libmetis_MUMPS_PrintTimers

/* util.c */
#define iargmax_strd                    libmetis_MUMPS_iargmax_strd 
#define iargmax_nrm                     libmetis_MUMPS_iargmax_nrm
#define iargmax2_nrm                    libmetis_MUMPS_iargmax2_nrm
#define rargmax2                        libmetis_MUMPS_rargmax2
#define InitRandom                      libmetis_MUMPS_InitRandom
#define metis_rcode                     libmetis_MUMPS_metis_rcode

/* wspace.c */
#define AllocateWorkSpace               libmetis_MUMPS_AllocateWorkSpace                  
#define AllocateRefinementWorkSpace     libmetis_MUMPS_AllocateRefinementWorkSpace
#define FreeWorkSpace                   libmetis_MUMPS_FreeWorkSpace
#define wspacemalloc                    libmetis_MUMPS_wspacemalloc
#define wspacepush                      libmetis_MUMPS_wspacepush
#define wspacepop                       libmetis_MUMPS_wspacepop
#define iwspacemalloc                   libmetis_MUMPS_iwspacemalloc
#define rwspacemalloc                   libmetis_MUMPS_rwspacemalloc
#define ikvwspacemalloc                 libmetis_MUMPS_ikvwspacemalloc
#define cnbrpoolReset                   libmetis_MUMPS_cnbrpoolReset
#define cnbrpoolGetNext                 libmetis_MUMPS_cnbrpoolGetNext
#define vnbrpoolReset                   libmetis_MUMPS_vnbrpoolReset
#define vnbrpoolGetNext                 libmetis_MUMPS_vnbrpoolGetNext

#endif


