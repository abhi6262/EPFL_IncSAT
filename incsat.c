/**CFile****************************************************************

  FileName    [incsat.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Iterative Layering - Incremental SAT Package.]

  Synopsis    [Checking Functional Dependency Using SAT solver.]

  Author      [Abhishek Sharma]
  
  Date        [Ver. 1.0. Started - May 17, 2013; updated July 1, 2013]

***********************************************************************/

#include "base/main/main.h"
#include "incsat.h"

ABC_NAMESPACE_IMPL_START

void genbricks_AbcNtk2AigNtkMapping(Abc_Ntk_t * pNtk, Cnf_Dat_t * pCnfA, char * pFileName);
Cnf_Dat_t * Cnf_DeriveSimple1( Aig_Man_t * p, int nOutputs, int * globnum, int * PiCnfNum );
Cnf_Dat_t * Cnf_DeriveSimple2( Abc_Ntk_t * pNtk, Aig_Man_t * p, int nOutputs, int * globnum, char ** PiName, int * PiCnfNum );
        
////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [Function for initializing the SAT problem]

  Description []
               
  SideEffects []

  SeeAlso     []

 ***********************************************************************/
/*       
sat_solver * IncSat_SatInit (Abc_Ntk_t * pNtk, Vec_Ptr_t * vBricks, int * globnum, char ** PiName, int * PiCnfNum, Abc_Ntk_t * pNtkMiterRes) {
    
    Abc_Ntk_t * pNtkMiterFunc;
    Abc_Ntk_t * pNtkMiterBricks;
    
    if (pNtk == NULL) {
        Abc_Print(-1, "The input network is empty.\n");
        return 0;        
    }
    
    if (!Abc_NtkIsStrash(pNtk)){
        Abc_Print(-1,"this command is only applicable to strashed networks.\n");
        return 0;
    }
    
    //Ana Petkovska's functions
    
    // Construct part of the network with the bricks
    pNtkMiterBricks = FD_CreateMiterBricks(vBricks);

    // Copy the network, make XOR of each output of the network and then connect them with a big OR gate
    pNtkMiterFunc = FD_CopyXorNetOutputsOr(pNtk);
    
    if(pNtkMiterBricks != NULL) {
        // Combine the miter from the network and the miter from the bricks
        pNtkMiterRes = FD_AndNet2ToNet1(pNtkMiterFunc, pNtkMiterBricks);
    } else {
        // if there were no bricks give to the SAT solver just the miter from the function
        pNtkMiterRes = pNtkMiterFunc;
    }
    
    Aig_Man_t * pAigMiter; 
    Cnf_Dat_t * pCnfA;
    Abc_Obj_t * pObj;
    int i, k;
    int * globnum;
    globnum = ABC_ALLOC(int , 1);
    
    pAigMiter = GenInt_NtkToDar(pNtkMiterRes);   //Ana's Function
    
    int * PiCnfNum; //required for adding brick to sat solver structure in the second function
    PiCnfNum = ABC_ALLOC( int, Aig_ManCiNum(pAigMiter) );
    
    char ** PiName; //required for adding brick to sat solver structure in the second function
    PiName = ABC_ALLOC( char *, Aig_ManCiNum(pAigMiter) );
    Abc_NtkForEachCi(pNtkMiterRes, pObj, k) {
            //strcpy(PiName[k], Abc_ObjName(pObj1));
        PiName[k] = Abc_ObjName(pObj);
        //Abc_Print(1, PiName[k]);
        //Abc_Print(1, "\n");
        }
    assert(Abc_NtkCiNum(pNtkMiterRes) == Aig_ManCiNum(pAigMiter) );
    
    
    pCnfA = Cnf_DeriveSimple1( pAigMiter, 1, globnum, PiCnfNum );
    genbricks_AbcNtk2AigNtkMapping(pNtkMiterRes, pCnfA, "print03_aigMapmiter.v");
    Cnf_DataWriteIntoFile( pCnfA, "mitercnf.txt", 0, NULL, NULL );
    
    char Buffer[500];
    sprintf(Buffer, "miter.v");
    Io_Write(pNtkMiterRes, Buffer, IO_FILE_VERILOG);
    
    sat_solver * pSat = NULL;
    pSat = sat_solver_new();
    sat_solver_setnvars(pSat, pCnfA->nVars ); 
    sat_solver_store_alloc(pSat);


    for (i = 0; i < pCnfA->nClauses; i++) {
        if (!sat_solver_addclause(pSat, pCnfA->pClauses[i], pCnfA->pClauses[i + 1])) {
            sat_solver_delete(pSat);
            //return NULL;
        }
    }
    
    sat_solver_store_write(pSat, "print06_satClausesoriginalmitercnf.txt");

    return pSat;
    
}





/**Function*************************************************************

  Synopsis    [Function for incrementing the SAT problem]

  Description []
               
  SideEffects []

  SeeAlso     []

 ***********************************************************************/
/*int IncSat_SatIncrement(sat_solver * satSolver, Abc_Ntk_t * pNtkBrick, int * globnum, char ** PiName, int * PiCnfNum) {
    
    Abc_Ntk_t * pNtkMiterBrick;
    pNtkMiterBrick = FD_CopyPlusXor(pNtkBrick, 1);   //Make two copies of the brick and XNor them
    
    Aig_Man_t * pAigMiterBrick; 
    Cnf_Dat_t * pCnfA;
    int i, flag;
    flag = 1;
    
    pAigMiterBrick = GenInt_NtkToDar(pNtkMiterBrick); 
    pCnfA = Cnf_DeriveSimple2( pNtkMiterBrick, pAigMiterBrick, 1, globnum, PiName, PiCnfNum );
    genbricks_AbcNtk2AigNtkMapping(pNtkMiterBrick, pCnfA, "print03_aigMapadditionalbrick.v");
    Cnf_DataWriteIntoFile( pCnfA, "additionalbrickcnf.txt", 0, NULL, NULL );
    
    
    for (i = 0; i < pCnfA->nClauses; i++) {
        if (!sat_solver_addclause(satSolver, pCnfA->pClauses[i], pCnfA->pClauses[i + 1])) {
            sat_solver_delete(satSolver);
            flag = 0;
            //return NULL;
        }
    }
    
    return flag;
    
    
}






/**Function*************************************************************

  Synopsis    [Function for solving the SAT problem]

  Description []
               
  SideEffects []

  SeeAlso     []

 ***********************************************************************/
/*Vec_Int_t * IncSat_SatSolve(sat_solver * satSolver, Vec_Int_t * vUsedBricks, Abc_Ntk_t * pNtkMiterRes) {
    
    Vec_Int_t * satassignment;
    Sto_Man_t * pCnfLearned;
    int RetValue, i;
    Abc_Obj_t * pObj;
    
    // solve the problem
    RetValue = sat_solver_solve( satSolver, NULL, NULL, (ABC_INT64_T)0, (ABC_INT64_T)0, (ABC_INT64_T)0, (ABC_INT64_T)0 );

    // report the result
    if ( RetValue == l_Undef )
    {
        Abc_Print(1, "The problem timed out.\n" );
    }
    
    else if ( RetValue == l_True )
    {
        Abc_Print(1, "The problem is SATISFIABLE.\n" );
        Abc_Print(1, "The target function does not functionally depend on the set of bricks.\n" );
    }
    
    else if ( RetValue == l_False )
    {
        Abc_Print(1, "The problem is UNSATISFIABLE.\n" );
        Abc_Print(1, "The target function functionally depends on the set of bricks.\n" );
    }
    
if ( RetValue == l_True )
    {
//        Vec_Int_t * vCiIds = Abc_NtkGetCiIds( pNtk );
        Vec_Int_t * vCiIds = Abc_NtkGetCiSatVarNums( pNtkMiterRes );
        pNtkMiterRes->pModel = Sat_SolverGetModel( satSolver, vCiIds->pArray, vCiIds->nSize );
        Vec_IntFree( vCiIds );
        FD_WriteCounter( pNtkMiterRes, pNtkMiterRes->pModel, 12 );
    }
    
    // get the learned clauses
    pCnfLearned = (Sto_Man_t *) sat_solver_store_release(satSolver);
    // delete the SAT solver
    sat_solver_delete(satSolver);
    
    if (pCnfLearned == NULL) {
        Abc_Print(-1, "Solving the SAT problem has failed.\n");
        return NULL;
    }
    
Sto_ManDumpClauses( pCnfLearned, "print06_satClauses12.txt" );

    satassignment = Vec_IntAlloc(Abc_NtkPiNum(pNtkMiterRes));
    Abc_NtkForEachPi( pNtkMiterRes, pObj, i )
        Vec_IntSetEntry(satassignment, i, pNtkMiterRes->pModel[i]);
    
    return satassignment;
    
 
}




/**Function*************************************************************

  Synopsis    [Print the mapping from Abc network to Aig network and varnum in the cnf]

  Description []
               
  SideEffects []

  SeeAlso     []

 ***********************************************************************/

/*void genbricks_AbcNtk2AigNtkMapping(Abc_Ntk_t * pNtk, Cnf_Dat_t * pCnfA, char * pFileName){
    Abc_Obj_t * pObj;
    int i;

    FILE * pFile = fopen(pFileName, "w");

    if (pFile == NULL) {
        fprintf(stdout, "AbcNtk2AigNtkMapping: Cannot open the output file \"%s\".\n", pFileName);
    } else {
        fprintf(pFile, "-- AbcNtk2AigNtkMapping: Written by ABC on %s\n", Extra_TimeStamp());
        fprintf(pFile, "\n");

        Abc_NtkForEachObj(pNtk, pObj, i) {
            fprintf(pFile, "objId: %d, objName: %s => aigNodeId %d (level: %d, type: %d", pObj->Id, Abc_ObjName(pObj), ((Aig_Obj_t *)pObj->pCopy)->Id, pObj->Level, pObj->Type);
            fprintf(pFile, ")");
            fprintf(pFile, " => cnf pVarNum: %d, literal: %d", pCnfA->pVarNums[((Aig_Obj_t *)pObj->pCopy)->Id], toLitCond(pCnfA->pVarNums[((Aig_Obj_t *)pObj->pCopy)->Id], 0));
            fprintf(pFile, "\n");
        }

        fprintf(pFile, "\n");
        fclose(pFile);
    }
}



/**Function*************************************************************

  Synopsis    [Cnf_DeriveSimple with globnum added for adding brick to a SAT solver structure.]

  Description [Also for remembering the varnums of inputs in an array.]
               
  SideEffects []

  SeeAlso     []

 ***********************************************************************/

/*Cnf_Dat_t * Cnf_DeriveSimple1( Aig_Man_t * p, int nOutputs, int * globnum, int * PiCnfNum )
{
    Aig_Obj_t * pObj;
    Cnf_Dat_t * pCnf;
    int OutVar, PoVar, pVars[32], * pLits, ** pClas;
    int i, nLiterals, nClauses, Number;

    // count the number of literals and clauses
    nLiterals = 1 + 7 * Aig_ManNodeNum(p) + Aig_ManCoNum( p ) + 3 * nOutputs;
    nClauses = 1 + 3 * Aig_ManNodeNum(p) + Aig_ManCoNum( p ) + nOutputs;

    // allocate CNF
    pCnf = ABC_ALLOC( Cnf_Dat_t, 1 );
    memset( pCnf, 0, sizeof(Cnf_Dat_t) );
    pCnf->pMan = p;
    pCnf->nLiterals = nLiterals;
    pCnf->nClauses = nClauses;
    pCnf->pClauses = ABC_ALLOC( int *, nClauses + 1 );
    pCnf->pClauses[0] = ABC_ALLOC( int, nLiterals );
    pCnf->pClauses[nClauses] = pCnf->pClauses[0] + nLiterals;

    // create room for variable numbers
    pCnf->pVarNums = ABC_ALLOC( int, Aig_ManObjNumMax(p) );
//    memset( pCnf->pVarNums, 0xff, sizeof(int) * Aig_ManObjNumMax(p) );
    for ( i = 0; i < Aig_ManObjNumMax(p); i++ )
        pCnf->pVarNums[i] = -1;
    // assign variables to the last (nOutputs) POs
    Number = 1;
    if ( nOutputs )
    {
//        assert( nOutputs == Aig_ManRegNum(p) );
//        Aig_ManForEachLiSeq( p, pObj, i )
//            pCnf->pVarNums[pObj->Id] = Number++;
        Aig_ManForEachCo( p, pObj, i )
            pCnf->pVarNums[pObj->Id] = Number++;
    }
    // assign variables to the internal nodes
    Aig_ManForEachNode( p, pObj, i )
        pCnf->pVarNums[pObj->Id] = Number++;
    
    // assign variables to the PIs and constant node
    Aig_ManForEachCi( p, pObj, i )
    {
        pCnf->pVarNums[pObj->Id] = Number++;
        PiCnfNum[i] = Number - 1;
    }
    pCnf->pVarNums[Aig_ManConst1(p)->Id] = Number++;
    pCnf->nVars = Number;
    (*globnum) = Number;
/*
    // print CNF numbers
    printf( "SAT numbers of each node:\n" );
    Aig_ManForEachObj( p, pObj, i )
        printf( "%d=%d ", pObj->Id, pCnf->pVarNums[pObj->Id] );
    printf( "\n" );
*/
    // assign the clauses
    /*pLits = pCnf->pClauses[0];
    pClas = pCnf->pClauses;
    Aig_ManForEachNode( p, pObj, i )
    {
        OutVar   = pCnf->pVarNums[ pObj->Id ];
        pVars[0] = pCnf->pVarNums[ Aig_ObjFanin0(pObj)->Id ];
        pVars[1] = pCnf->pVarNums[ Aig_ObjFanin1(pObj)->Id ];

        // positive phase
        *pClas++ = pLits;
        *pLits++ = 2 * OutVar; 
        *pLits++ = 2 * pVars[0] + !Aig_ObjFaninC0(pObj); 
        *pLits++ = 2 * pVars[1] + !Aig_ObjFaninC1(pObj); 
        // negative phase
        *pClas++ = pLits;
        *pLits++ = 2 * OutVar + 1; 
        *pLits++ = 2 * pVars[0] + Aig_ObjFaninC0(pObj); 
        *pClas++ = pLits;
        *pLits++ = 2 * OutVar + 1; 
        *pLits++ = 2 * pVars[1] + Aig_ObjFaninC1(pObj); 
    }
 
    // write the constant literal
    OutVar = pCnf->pVarNums[ Aig_ManConst1(p)->Id ];
    assert( OutVar <= Aig_ManObjNumMax(p) );
    *pClas++ = pLits;
    *pLits++ = 2 * OutVar;  

    // write the output literals
    Aig_ManForEachCo( p, pObj, i )
    {
        OutVar = pCnf->pVarNums[ Aig_ObjFanin0(pObj)->Id ];
        if ( i < Aig_ManCoNum(p) - nOutputs )
        {
            *pClas++ = pLits;
            *pLits++ = 2 * OutVar + Aig_ObjFaninC0(pObj); 
        }
        else
        {
            PoVar  = pCnf->pVarNums[ pObj->Id ];
            // first clause
            *pClas++ = pLits;
            *pLits++ = 2 * PoVar; 
            *pLits++ = 2 * OutVar + !Aig_ObjFaninC0(pObj); 
            // second clause
            *pClas++ = pLits;
            *pLits++ = 2 * PoVar + 1; 
            *pLits++ = 2 * OutVar + Aig_ObjFaninC0(pObj); 
        }
    }

    // verify that the correct number of literals and clauses was written
    assert( pLits - pCnf->pClauses[0] == nLiterals );
    assert( pClas - pCnf->pClauses == nClauses );
    return pCnf;
}




/**Function*************************************************************

  Synopsis    [Cnf_DeriveSimple with globnum to derive the cnf of a brick to be added to an existing SAT solver structure.]

  Description []
               
  SideEffects []

  SeeAlso     []

 ***********************************************************************/

/*Cnf_Dat_t * Cnf_DeriveSimple2( Abc_Ntk_t * pNtk, Aig_Man_t * p, int nOutputs, int * globnum, char ** PiName, int * PiCnfNum )
{
    Aig_Obj_t * pObj;
    Cnf_Dat_t * pCnf;
    int OutVar, PoVar, pVars[32], * pLits, ** pClas;
    int i, k, nLiterals, nClauses, Number, Number1;

    // count the number of literals and clauses
    nLiterals = 1 + 7 * Aig_ManNodeNum(p) + Aig_ManCoNum( p ) + 3 * nOutputs;
    nClauses = 1 + 3 * Aig_ManNodeNum(p) + Aig_ManCoNum( p ) + nOutputs;

    // allocate CNF
    pCnf = ABC_ALLOC( Cnf_Dat_t, 1 );
    memset( pCnf, 0, sizeof(Cnf_Dat_t) );
    pCnf->pMan = p;
    pCnf->nLiterals = nLiterals;
    pCnf->nClauses = nClauses;
    pCnf->pClauses = ABC_ALLOC( int *, nClauses + 1 );
    pCnf->pClauses[0] = ABC_ALLOC( int, nLiterals );
    pCnf->pClauses[nClauses] = pCnf->pClauses[0] + nLiterals;

    // create room for variable numbers
    pCnf->pVarNums = ABC_ALLOC( int, Aig_ManObjNumMax(p) );
//    memset( pCnf->pVarNums, 0xff, sizeof(int) * Aig_ManObjNumMax(p) );
    for ( i = 0; i < Aig_ManObjNumMax(p); i++ )
        pCnf->pVarNums[i] = -1;
    // assign variables to the last (nOutputs) POs
    Number = (*globnum);
    Number1 = 1;
    if ( nOutputs )
    {
//        assert( nOutputs == Aig_ManRegNum(p) );
//        Aig_ManForEachLiSeq( p, pObj, i )
//            pCnf->pVarNums[pObj->Id] = Number++;
        Aig_ManForEachCo( p, pObj, i )
        { pCnf->pVarNums[pObj->Id] = Number++; Number1++; }
    }
    // assign variables to the internal nodes
    Aig_ManForEachNode( p, pObj, i )
    { pCnf->pVarNums[pObj->Id] = Number++; Number1++; }
    // assign variables to the PIs and constant node
    Aig_ManForEachCi( p, pObj, i )
    {
        Number1++;
        for ( k=0; k < (Abc_NtkCiNum(pNtk)); k++)
        {
            if ( strcmp(Abc_ObjName(Abc_NtkCi( pNtk, i )),PiName[k]) == 0 ) {pCnf->pVarNums[pObj->Id] = PiCnfNum[k];}
            //Abc_Print(1, Abc_ObjName(Abc_NtkCi( pNtk, i ))); Abc_Print(1, PiName[k]); Abc_Print(1, PiCnfNum[k]); Abc_Print(1, "\n");}
        }
        //Abc_Print(1, "\n");
    }
    pCnf->pVarNums[Aig_ManConst1(p)->Id] = Number++;
    pCnf->nVars = ++Number1;
    (*globnum) = Number;
/*
    // print CNF numbers
    printf( "SAT numbers of each node:\n" );
    Aig_ManForEachObj( p, pObj, i )
        printf( "%d=%d ", pObj->Id, pCnf->pVarNums[pObj->Id] );
    printf( "\n" );
*/
    // assign the clauses
   /* pLits = pCnf->pClauses[0];
    pClas = pCnf->pClauses;
    Aig_ManForEachNode( p, pObj, i )
    {
        OutVar   = pCnf->pVarNums[ pObj->Id ];
        pVars[0] = pCnf->pVarNums[ Aig_ObjFanin0(pObj)->Id ];
        pVars[1] = pCnf->pVarNums[ Aig_ObjFanin1(pObj)->Id ];

        // positive phase
        *pClas++ = pLits;
        *pLits++ = 2 * OutVar; 
        *pLits++ = 2 * pVars[0] + !Aig_ObjFaninC0(pObj); 
        *pLits++ = 2 * pVars[1] + !Aig_ObjFaninC1(pObj); 
        // negative phase
        *pClas++ = pLits;
        *pLits++ = 2 * OutVar + 1; 
        *pLits++ = 2 * pVars[0] + Aig_ObjFaninC0(pObj); 
        *pClas++ = pLits;
        *pLits++ = 2 * OutVar + 1; 
        *pLits++ = 2 * pVars[1] + Aig_ObjFaninC1(pObj); 
    }
 
    // write the constant literal
    OutVar = pCnf->pVarNums[ Aig_ManConst1(p)->Id ];
    //assert( OutVar <= Aig_ManObjNumMax(p) );
    *pClas++ = pLits;
    *pLits++ = 2 * OutVar;  

    // write the output literals
    Aig_ManForEachCo( p, pObj, i )
    {
        OutVar = pCnf->pVarNums[ Aig_ObjFanin0(pObj)->Id ];
        if ( i < Aig_ManCoNum(p) - nOutputs )
        {
            *pClas++ = pLits;
            *pLits++ = 2 * OutVar + Aig_ObjFaninC0(pObj); 
        }
        else
        {
            PoVar  = pCnf->pVarNums[ pObj->Id ];
            // first clause
            *pClas++ = pLits;
            *pLits++ = 2 * PoVar; 
            *pLits++ = 2 * OutVar + !Aig_ObjFaninC0(pObj); 
            // second clause
            *pClas++ = pLits;
            *pLits++ = 2 * PoVar + 1; 
            *pLits++ = 2 * OutVar + Aig_ObjFaninC0(pObj); 
        }
    }

    // verify that the correct number of literals and clauses was written
    assert( pLits - pCnf->pClauses[0] == nLiterals );
    assert( pClas - pCnf->pClauses == nClauses );
    return pCnf;
}

*/


        
////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END
