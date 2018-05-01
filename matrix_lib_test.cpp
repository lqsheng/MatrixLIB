#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include "Matrix.h"

using std::cout;
using std::endl;

//对于较大的矩阵,用此函数打印它的缩写
template  <class ElemType>
void print_large_matrix(Matrix<ElemType> A)
{
    if((A.getCols()>1)&&(A.getRows()>1)){
        std::cout<<std::setw(10)<<A(0,0)<<" , ... , "<<std::setw(10)<<A(0,A.getCols()-1)<<std::endl;
        std::cout<<"    :      , ... ,     :      "<<"size="<<A.getRows()<<"x"<<A.getCols()<<std::endl;
        std::cout<<std::setw(10)<<A(A.getRows()-1,0)<<" , ... , "<<std::setw(10)<<A(A.getRows()-1,A.getCols()-1)<<std::endl;
    }
    else if(A.getRows()>1){
        std::cout<<std::setw(10)<<A(0,0)<<std::endl;
        std::cout<<"    :      "<<"size="<<A.getRows()<<"x"<<A.getCols()<<std::endl;
        std::cout<<std::setw(10)<<A(A.getRows()-1,0)<<std::endl;
    }
    else if(A.getCols()>1){
        std::cout<<std::setw(10)<<A(0,0)<<" , ... , "<<std::setw(10)<<A(0,A.getCols()-1)<<"size="<<A.getRows()<<"x"<<A.getCols()<<std::endl;
    }
}

int main(int argc, char *argv[])
{
    Matrix<int> m_a(5,5,1);
    Matrix<int> m_b(5,2);
    Matrix<int> m_c=m_b-m_a;
    cout<<"m_a="<<endl<<m_a;
    cout<<"m_b="<<endl<<m_b;
    cout<<"m_b-m_a="<<endl<<m_c;
    cout<<"m_b+m_a="<<endl<<m_b+m_a;
    cout<<"m_b+m_a+m_a="<<endl<<m_b+m_a+m_a;
    cout<<"m_b*m_a="<<endl<<m_b*m_a;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Matrix<double> m_d(100,100,0);
    Matrix<double> m_e(100,100,0);
    for(int i=0;i<100;i++){//生成随机矩阵
        for(int j=0;j<100;j++){
                m_d(i,j)=u(e);
        }
    }
    cout<<"m_d="<<endl;
    print_large_matrix(m_d);
    m_e=m_d.transpose()*m_d;//构造正定的随机A矩阵
    cout<<"m_e="<<endl;
    print_large_matrix(m_e);
    cout<<"start to get inverse matrix"<<endl;
    clock_t time_start=clock();
    Matrix<double> m_f=m_e.inverse();
    std::cout<<"time used is "<<1000*(clock()-time_start)/(double)CLOCKS_PER_SEC<<"ms"<<std::endl;
    cout<<"m_e.inverse()="<<endl;
    print_large_matrix(m_f);
    cout<<"m_e*m_e.inverse()="<<endl;
    print_large_matrix(m_e*m_f);
    return 0;
}
