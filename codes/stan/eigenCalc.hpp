// namespace model_SEIR_namespace {

template <typename T>
inline 
    // typename boost::math::tools::promote_args<T>::type
    typename boost::math::tools::promote_args<T>::type
    maxEigen(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M, std::ostream* pstream__){
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > solver(M, Eigen::EigenvaluesOnly);
    // int n = M.rows();/
    Eigen::Matrix<T, Eigen::Dynamic, 1> e = solver.eigenvalues().real();
    return e.maxCoeff();;
    // return .2;
}
// }