#ifndef libmatrix
#define libmatrix
#ifdef WITH_NEON
#include <arm_neon.h>
#endif
#include </opt/homebrew/opt/libomp/include/omp.h>

//智能指针类
template<typename T>
class SmartPtr{//智能指针类
    struct MatData{
        size_t row;//完整2D矩阵大小
        size_t col;
        T* data;
        int cnt;//引用次数
        ~MatData(){
            try {
                delete [] data;
                data=nullptr;
            }catch(...){}
        }
    };

    public:
    MatData * totalData;
    //构造器
    SmartPtr(size_t row,size_t col):totalData(nullptr){
        try {
            totalData=new MatData{row,col,new T[row*col],1};
        }catch(...){
            try{delete[] totalData->data; totalData->data=nullptr;} catch(...){}
            try{delete totalData;totalData=nullptr;}catch(...){}
            throw;
        }
    }
    //析构函数
    ~SmartPtr(){
        try{
            if(--totalData->cnt==0){
                delete totalData;
                totalData=nullptr;
            }
        }catch(...){}
    }
    //拷贝构造器
    SmartPtr(const SmartPtr& p){
        this->totalData=p.totalData;
        ++this->totalData->cnt;
    }
    //重载赋值运算符
    SmartPtr &operator=(const SmartPtr &p){
        if(this==&p || this->totalData==p.totalData) return *this;
        //删掉自己原来的数据
        if(--totalData->cnt==0) delete totalData;
        totalData = p.totalData;
        ++totalData->cnt;
        return *this;
    }
    //取任意位置的值
    T &operator()(size_t r,size_t c){
        return totalData->data[r * totalData->col + c];
    }
};

template <typename T>
class Matrix{
    public:
    size_t rows,cols,drow,dcol;
    std::vector<SmartPtr<T>> data;
    //构造器
    Matrix(size_t r,size_t c,int ch=1):rows(r),cols(c),drow(0),dcol(0){
        for(int i=0;i<ch;i++){
            data.emplace_back(r,c);
        }
    }
    //设置ROI大小
    Matrix resize(size_t r,size_t c,size_t dr=0,size_t dc=0){
        try{
        if(dr+r>this->data.at(0).totalData->row) throw("row index out of range");
        if(dc+c>this->data.at(0).totalData->col) throw("column index out of range");
        this->rows=r;
        this->cols=c;
        this->drow=dr;
        this->dcol=dc;
        return *this;}
        catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
    }
    Matrix<T> operator+(Matrix a) {return addMatrix(*this,a);}
    Matrix<T> operator*(Matrix a) {return multiplyMatrix(*this,a);}
    Matrix<T> operator-(Matrix a) {return subtractMatrix(*this,a);}
    Matrix<T> operator*(T n) {return multiplyMatrix(*this,n);}
    Matrix<T> operator+(T n) {return addMatrix(*this,n);}
    Matrix<T> operator-(T n) {return subtractMatrix(*this,n);}
};


//输入矩阵
template <typename T>
void fillMatrix(Matrix<T> a){
    size_t rows=a.rows,cols=a.cols;
    int channels=a.data.size();
    size_t totalcols=a.data.at(0).totalData->col;
    for(int i=0;i<channels;i++){
        printf("Please input a %zu by %zu matrix\n",rows,cols);
        for(size_t j=0;j<rows;j++){
            for(size_t k=0;k<cols;k++){
                cin>>a.data.at(i).totalData->data[(j+a.drow)*totalcols+k+a.dcol];
            }
        }
    }
    cout<<"The input matrix is:"<<endl;
    printMatrix(a);
    cout<<endl;
};

//打印矩阵
template <typename T>
void printMatrix(const Matrix<T> a){
    try{
        if(a.data.size()==0 ) throw("This matrix is empty");
        int channels=a.data.size();
        size_t totalcols=a.data.at(0).totalData->col;
        for(int c=0;c<channels;c++){
            if(a.data.at(c).totalData->data==nullptr) throw("This matrix is empty");
            printf("Channel %d:\n",c);
            for(size_t i=a.drow;i<a.drow+a.rows;i++){
                printf("[ ");
                for(size_t j=a.dcol;j<a.dcol+a.cols;j++){
                    cout<<a.data.at(c).totalData->data[i*totalcols+j]<<" ";
                }
                printf("]\n");
            }
        }
    }
    catch(const char* msg) {
        cerr<<msg<<endl;
    }
}

//深拷贝矩阵
template <typename T>
Matrix<T> copyMatrix(const Matrix<T> a){
    try{
    if(a.data.size()==0 ) throw("This matrix is empty");
    Matrix<T> result(a.rows,a.cols,a.data.size());
    int channels=a.data.size();
    size_t totalcols=a.data.at(0).totalData->col;
    for(int c=0;c<channels;c++){
        if(a.data.at(c).totalData->data==nullptr) throw("This matrix is empty");
        for(size_t i=0;i<a.rows;i++){
            for(size_t j=0;j<a.cols;j++){
                result.data.at(c).totalData->data[i*a.cols+j]=a.data.at(c).totalData->data[(i+a.drow)*totalcols+j+a.dcol];
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
}
//生成随机矩阵
template <typename T>
Matrix<T> createMatrixRandom(int channels,const size_t rows, const size_t cols, const T leftBound, const T rightBound) {
    Matrix<T> newMatrix(rows, cols,channels);
    for(int c=0;c<channels;c++){
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++) {
                newMatrix.data[i * cols + j] =
                leftBound + (T)((rightBound - leftBound) * (float) rand() / (RAND_MAX + 1.0));
            }
        }
    }
    return newMatrix;
}
//矩阵加矩阵
template <typename T>
Matrix<T> addMatrix(const Matrix<T> a,const Matrix<T> b){
    try{
    if(a.data.size()==0 ) throw("This matrix is empty");
    if(a.data.size()!=b.data.size()) throw("#Channels of mat1 and mat2 doesn't match");
    Matrix result=copyMatrix(a);
    int channels=a.data.size();
    size_t bcols=b.data.at(0).totalData->col;
    for(int c=0;c<channels;c++){
        if(b.data.at(c).totalData->data==nullptr) throw("This matrix is empty");
        for(size_t i=0;i<a.rows;i++){
            for(size_t j=0;j<a.cols;j++){
                result.data.at(c).totalData->data[i*a.cols+j]+=b.data.at(c).totalData->data[(i+b.drow)*bcols+j+b.dcol];
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
};
//矩阵加上一个数
template <typename T>
Matrix<T> addMatrix(const Matrix<T> a,T n){
    try{
    if(a.data.size()==0 ) throw("This matrix is empty");
    Matrix result=copyMatrix(a);
    int channels=a.data.size();
    size_t totalcols=a.data.at(0).totalData->col;
    for(int c=0;c<channels;c++){
        for(size_t i=0;i<a.rows;i++){
            for(size_t j=0;j<a.cols;j++){
                result.data.at(c).totalData->data[i*a.cols+j]+=n;
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
};
//矩阵减一个数
template <typename T>
Matrix<T> subtractMatrix(const Matrix<T> a,T n){
    try{
    if(a.data.size()==0 ) throw("This matrix is empty");
    Matrix result=copyMatrix(a);
    int channels=a.data.size();
    size_t totalcols=a.data.at(0).totalData->col;
    for(int c=0;c<channels;c++){
        for(size_t i=0;i<a.rows;i++){
            for(size_t j=0;j<a.cols;j++){
                result.data.at(c).totalData->data[i*a.cols+j]-=n;
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
};
//矩阵减矩阵
template <typename T>
Matrix<T> subtractMatrix(const Matrix<T> a,const Matrix<T> b){
    try{
    if(a.data.size()==0 ) throw("This matrix is empty");
    if(a.data.size()!=b.data.size()) throw("#Channels of mat1 and mat2 doesn't match");
    Matrix result=copyMatrix(a);
    int channels=a.data.size();
    size_t bcols=b.data.at(0).totalData->col;
    for(int c=0;c<channels;c++){
        if(b.data.at(c).totalData->data==nullptr) throw("This matrix is empty");
        for(size_t i=0;i<a.rows;i++){
            for(size_t j=0;j<a.cols;j++){
                result.data.at(c).totalData->data[i*a.cols+j]-=b.data.at(c).totalData->data[(i+b.drow)*bcols+j+b.dcol];
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
};
//矩阵乘以一个数
template <typename T>
Matrix<T> multiplyMatrix(const Matrix<T> a,T n){
    try{
    if(a.data.size()==0 ) throw("This matrix is empty");
    Matrix result=copyMatrix(a);
    int channels=a.data.size();
    size_t totalcols=a.data.at(0).totalData->col;
    for(int c=0;c<channels;c++){
        for(size_t i=0;i<a.rows;i++){
            for(size_t j=0;j<a.cols;j++){
                result.data.at(c).totalData->data[i*a.cols+j]*=n;
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
};
//矩阵转置
template <typename T>
Matrix<T> transposeMatrix(const Matrix<T> a){
    try{
    if(a.data.size()==0 ) throw("This matrix is empty");
    size_t rows=a.rows,cols=a.cols;
    int channels=a.data.size();
    size_t totalcols=a.data.at(0).totalData->col;
    Matrix<T> result(cols,rows,channels);
    for(int ch=0;ch<channels;ch++){
        if(a.data.at(ch).totalData->data==nullptr) throw("This matrix is empty");
        T* p1=&a.data.at(ch).totalData->data[a.dcol+totalcols*a.drow];
        T* p2=&result.data.at(ch).totalData->data[0];
        for(size_t i=0;i<rows;i++){
            p1=&a.data.at(ch).totalData->data[a.dcol+totalcols*(i+a.drow)];
            p2=&result.data.at(ch).totalData->data[i];
            for(size_t j=0;j<cols;j++){
                *(p2)=*(p1++);
                p2+=rows;
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
}

//矩阵乘矩阵
template <typename T>
Matrix<T> multiplyMatrix(const Matrix<T> a,const Matrix<T> b){
    try{
    if(a.data.size()==0 ||b.data.size()==0) throw("This matrix is empty");
    size_t rows=a.rows;
    size_t cols=b.cols;
    size_t common=a.cols;
    if(common!=b.rows) throw("#Rows of mat1 doesn't match #cols of mat2");
    size_t acol=a.data.at(0).totalData->col;
    size_t bcol=b.data.at(0).totalData->col;
    int channels=a.data.size();
    if(channels!=b.data.size()) throw("#Channels of mat1 and mat2 doesn't match");
    Matrix<T> result(rows,cols,channels);
    #pragma omp parallel for
    for(int ch=0;ch<channels;ch++){
        for(size_t m=0;m<rows;m++){
            for(size_t n=0;n<cols;n++){
                T cache=0;
                for (size_t i=0;i<common;i++){
                    cache+=a.data.at(ch).totalData->data[(m+a.drow)*acol+i+a.dcol]*b.data.at(ch).totalData->data[(i+b.drow)*bcol+n+b.dcol];
                }
                result.data.at(ch).totalData->data[m*cols+n]=cache;
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
};

//矩阵乘矩阵（优化）
template <>
Matrix<float> multiplyMatrix(const Matrix<float> mata,const Matrix<float> matb){
    try{
    if(mata.data.size()==0 ||matb.data.size()==0) throw("This matrix is empty");
    Matrix<float> mat1=copyMatrix(mata);
    size_t channels=mat1.data.size();
    size_t rows=mat1.rows;
    size_t cols=matb.cols;
    size_t common=mat1.cols;
    if(common!=matb.rows) throw("#Rows of mat1 doesn't match #cols of mat2");
    Matrix<float> result(rows,cols,channels);
    size_t common_mod4=common-(common%4);
    Matrix<float> mat2=transposeMatrix(matb);
    #pragma omp parallel for
    for(size_t ch=0;ch<channels;ch++){
        for(size_t m=0;m<rows;m++){
            float * p3=&(result.data.at(ch).totalData->data[m*cols]);
            for(size_t n=0;n<cols;n++){
                float * p1=&mat1.data.at(ch).totalData->data[m * common];
                float * p2=&mat2.data.at(ch).totalData->data[n * common];
                float32x4_t a, b;
                float32x4_t c = vdupq_n_f32(0);
                float sum[4] = {0};
                //是4倍数的部分
                for (size_t i=0;i<common_mod4;i+=4){
                    a = vld1q_f32(p1);
                    b = vld1q_f32(p2);
                    c = vaddq_f32(c, vmulq_f32(a, b));
                    p1 += 4;
                    p2 += 4;
                }
                vst1q_f32(sum, c);
                *p3=(sum[0]+sum[1]+sum[2]+sum[3]);
                //除4剩余的部分
                for (size_t i=common_mod4;i<common;i++){
                    *p3 += *(p1++) * *(p2++);
                }
                p3++;
            }
        }
    }
    return result;}
    catch(const char* msg) {
        cerr<<msg<<endl;
        throw(msg);
    }
};


#endif
