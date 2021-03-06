#include "TpetraUtil.hpp"

#include <limits>

void dumpTCMAT(const Teuchos::RCP<const TCMAT> &A, std::string filename) {
    filename = filename + std::string("_TCMAT.mtx");
    if (A->getComm()->getRank() == 0) {
        std::cout << "dumping" << filename << std::endl;
        A->print(std::cout);
    }

    Tpetra::MatrixMarket::Writer<TCMAT> matDumper;
    matDumper.writeSparseFile(filename, A, filename, filename, true);
}

void dumpTMV(const Teuchos::RCP<const TMV> &A, std::string filename) {
    filename = filename + std::string("_TMV.mtx");
    if (A->getMap()->getComm()->getRank() == 0) {
        std::cout << "dumping" << filename << std::endl;
        A->print(std::cout);
    }

    Tpetra::MatrixMarket::Writer<TMV> matDumper;
    matDumper.writeDenseFile(filename, A, filename, filename);
}

void dumpTV(const Teuchos::RCP<const TV> &A, std::string filename) {
    filename = filename + std::string("_TV.mtx");
    if (A->getMap()->getComm()->getRank() == 0) {
        std::cout << "dumping" << filename << std::endl;
        A->print(std::cout);
    }

    Tpetra::MatrixMarket::Writer<TV> matDumper;
    matDumper.writeDenseFile(filename, A, filename, filename);
}

void dumpTOP(const Teuchos::RCP<const TOP> &A, std::string filename) {
    filename = filename + std::string("_TOP.mtx");

    auto domainMapRcp = A->getDomainMap();
    auto rangeMapRcp = A->getRangeMap();
    if (domainMapRcp->getComm()->getRank() == 0) {
        std::cout << "dumping" << filename << std::endl;
        std::cout << A->description() << std::endl;
    }
    Teuchos::RCP<TMV> eyeRcp = Teuchos::rcp(new TMV(domainMapRcp, domainMapRcp->getGlobalNumElements(), true));
    Teuchos::RCP<TMV> resultRcp = Teuchos::rcp(new TMV(rangeMapRcp, domainMapRcp->getGlobalNumElements(), true));

    // fill eyeRcp as identity
    for (int c = 0; c < eyeRcp->getNumVectors(); c++) {
        eyeRcp->replaceGlobalValue(c, c, 1.0);
    }

    domainMapRcp->getComm()->barrier();
    A->apply(*eyeRcp, *resultRcp);

    Tpetra::MatrixMarket::Writer<TMV> matDumper;
    matDumper.writeDenseFile(filename, resultRcp, filename, filename);
}

Teuchos::RCP<const TCOMM> getMPIWORLDTCOMM() { return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD)); }

// return a fully copied TMAP with a given global size
Teuchos::RCP<TMAP> getFullCopyTMAPFromGlobalSize(const int &globalSize, Teuchos::RCP<const TCOMM> &commRcp) {
    std::vector<int> globalIndexOnLocal(globalSize);
#pragma omp parallel for schedule(dynamic, 2048)
    // define the matrix with col map = full cols on every node
    for (int kk = 0; kk < globalSize; kk++) {
        globalIndexOnLocal[kk] = kk;
    }
    return Teuchos::rcp(new TMAP(globalSize, globalIndexOnLocal.data(), globalSize, 0, commRcp));
}

// return a contiguous TMAP from local Size
Teuchos::RCP<TMAP> getTMAPFromLocalSize(const int &localSize, Teuchos::RCP<const TCOMM> &commRcp) {
    int globalSize = localSize;
    MPI_Allreduce(MPI_IN_PLACE, &globalSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return Teuchos::rcp(new TMAP(globalSize, localSize, 0, commRcp));
}

Teuchos::RCP<TV> getTVFromVector(const std::vector<double> &in, Teuchos::RCP<const TCOMM> &commRcp) {

    const int localSize = in.size();

    Teuchos::RCP<TMAP> contigMapRcp = getTMAPFromLocalSize(localSize, commRcp);

    Teuchos::RCP<TV> out = Teuchos::rcp(new TV(contigMapRcp, false));

    auto out_2d = out->getLocalView<Kokkos::HostSpace>();
    assert(out_2d.dimension_0() == localSize);

    out->modify<Kokkos::HostSpace>();
    for (int c = 0; c < out_2d.dimension_1(); c++) {
#pragma omp parallel for schedule(dynamic, 1024)
        for (int i = 0; i < out_2d.dimension_0(); i++) {
            out_2d(i, c) = in[i];
        }
    }

    return out;
}

Teuchos::RCP<TMV> getTMVFromVector(const std::vector<std::vector<double>> &in, Teuchos::RCP<const TCOMM> &commRcp) {
    const int nCol = in.size();
    // look for the minimal local size of in vecs
    int localSize = std::numeric_limits<int>::max();
    for (const auto &v : in) {
        localSize = std::min(localSize, (int)v.size());
    }

    Teuchos::RCP<TMAP> contigMapRcp = getTMAPFromLocalSize(localSize, commRcp);

    Teuchos::RCP<TMV> out = Teuchos::rcp(new TMV(contigMapRcp, nCol, false));

    auto out_2d = out->getLocalView<Kokkos::HostSpace>();
    assert(out_2d.dimension_0() == localSize);
    assert(out_2d.dimension_1() == nCol);

    out->modify<Kokkos::HostSpace>();
    for (int c = 0; c < out_2d.dimension_1(); c++) {
#pragma omp parallel for schedule(dynamic, 1024)
        for (int i = 0; i < out_2d.dimension_0(); i++) {
            out_2d(i, c) = in[c][i];
        }
    }
    return out;
}
