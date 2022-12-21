
#include <fstream>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static
void
PrintUsage (const char* progName)
{
    Print() << std::endl
            << "This utility performs a diff operation between specific"      << std::endl
            << "components of two"           << std::endl
            << "plotfiles that have the same geometrical domain and nodality" << std::endl
            << "(supports all nodality types; cell, face, edge, node)"        << std::endl
            << "and outputs the L0, L1, and L2 norms"                         << std::endl
            << "L1 = sum(|diff_ijk|)/npts"                       << std::endl
            << "L2 = sqrt[sum(diff_ijk^2)]/sqrt(npts)"           << std::endl
            <<  "(only single-level supported)"                               << std::endl << std::endl;

    Print() << "Usage:" << '\n';
    Print() << progName << '\n';
    Print() << "    infile1 = inputFileName1" << '\n';
    Print() << "      comp1 = component1    " << '\n';
    Print() << "    infile2 = inputFileName2" << '\n';
    Print() << "      comp2 = component2    " << '\n';        
    Print() << "    diffile = differenceFileName" << '\n';
    Print() << "              (If not specified no file is written)" << '\n' << '\n';

    Print() << "You can either point to the plotfile base directory itself, e.g."      << std::endl
            << "  infile=plt00000"                                                     << std::endl
            << "Or the raw data itself, e.g."                                          << std::endl
            << "  infile=plt00000/Level_0/Cell"                                        << std::endl
            << "the latter is useful for some applications that dump out raw"          << std::endl
            << "nodal data within a plotfile directory."                               << std::endl
            << "The program will first try appending 'Level_0/Cell'"                   << std::endl
            << "onto the specified filenames."                                         << std::endl
            << "If that _H file doesn't exist, it tries using the full specified name" << std::endl << std::endl;

    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {

        if (argc == 1) {
            PrintUsage(argv[0]);
        }

        const std::string farg = amrex::get_command_argument(1);
        if (farg == "-h" || farg == "--help")
        {
            PrintUsage(argv[0]);
        }

        // plotfile names for the coarse, fine, and subtracted output
        std::string iFile1, iFile2, difFile="";

        // read in parameters from inputs file
        ParmParse pp;

        // coarse MultiFab
        pp.query("infile1", iFile1);
        if (iFile1.empty())
            amrex::Abort("You must specify `infile1'");

        // fine MultiFab (might have same resolution as coarse)
        pp.query("infile2", iFile2);
        if (iFile2.empty())
            amrex::Abort("You must specify `infile2'");

        // subtracted output (optional)
        pp.query("diffile", difFile);

        int comp1=0;
        int comp2=0;
        pp.query("comp1",comp1);
        pp.query("comp2",comp2);

        // single-level for now
        // AMR comes later, where we iterate over each level in isolation

        // check to see whether the user pointed to the plotfile base directory
        // or the data itself
        if (amrex::FileExists(iFile1+"/Level_0/Cell_H")) {
            iFile1 += "/Level_0/Cell";
        }
        if (amrex::FileExists(iFile2+"/Level_0/Cell_H")) {
            iFile2 += "/Level_0/Cell";
        }

        // storage for the input coarse and fine MultiFabs
        MultiFab mf_1, mf_2;

        // read in plotfiles, 'coarse' and 'fine' to MultiFabs
        // note: fine could be the same resolution as coarse
        VisMF::Read(mf_1, iFile1);
        VisMF::Read(mf_2, iFile2);

        if (mf_1.contains_nan()) {
            Abort("First plotfile contains NaN(s)");
        }
        if (mf_2.contains_nan()) {
            Abort("Second plotfile contains NaN(s)");
        }

        // check nodality
        IntVect c_nodality = mf_1.ixType().toIntVect();
        IntVect f_nodality = mf_2.ixType().toIntVect();
        if (c_nodality != f_nodality) {
            Abort("plotfiles do not have the same nodality");
        }

        // get BoxArray and DistributionMapping
        BoxArray ba = mf_1.boxArray();
        DistributionMapping dm = mf_1.DistributionMap();

        // minimalBox() computes a single box to enclose all the boxes
        // enclosedCells() converts it to a cell-centered Box
        Box bx = ba.minimalBox().enclosedCells();
 
        // number of cells in the domain
        Print() << "npts in domain = " << bx.numPts() << std::endl;
        long npts = bx.numPts();
        
        // subtract coarse from coarsened fine
        MultiFab::Subtract(mf_1,mf_2,comp2,comp1,1,0);

        // force periodicity so faces/edges/nodes get weighted accordingly for L1 and L2 norms
        IntVect iv(AMREX_D_DECL(bx.length(0),
                                bx.length(1),
                                bx.length(2)));
        Periodicity period(iv);

        // compute norms of difference
        Real norm0 = mf_1.norm0(comp1);
        Real norm1 = mf_1.norm1(comp1,period);
        Real norm2 = mf_1.norm2(comp1,period);
        Print() << "(L0,L1,L2) " << norm0 << " " << norm1/npts << " " << norm2/sqrt(npts) << " " << std::endl;

        // write out the subtracted plotfile if diffile was specified at the command line
        if (difFile != "") {

            // define the problem domain as (0,1) for now
            RealBox real_box({AMREX_D_DECL(0.,0.,0.)},
                             {AMREX_D_DECL(1.,1.,1.)});

            Vector<int> is_periodic(AMREX_SPACEDIM,1);

            // build a geometry object so we can use WriteSingleLevelPlotfile
            Geometry geom(bx,&real_box,CoordSys::cartesian,is_periodic.data());

            MultiFab diff(ba,dm,1,0);
            MultiFab::Copy(diff,mf_1,0,comp1,1,0);            
            
            WriteSingleLevelPlotfile(difFile,diff,{"diff"},geom,0.,0);
        }

    }
    amrex::Finalize();
}
