#include <iostream>
#include <cstdlib> 
#include <vector> 
#include <sstream>
#include <fstream>
#include <chrono> 
#include <exception>
#include <math.h>

using uint = unsigned int;
using std::chrono::time_point;

struct Pair
{
    long loc1;
    long loc2;
    double dist;
};

class ConsoleCls {
public:
    template<typename T> 
    void Write(const T& out) {
        std::cout << out;
    }

    template<typename T> 
    void WriteLine(const T& out) {
        std::cout << out << std::endl;
    }

    void WriteLine() {
        std::cout << std::endl;
    }

    void WriteLine(const char c) {
        std::cout << std::string(&c) << std::endl;
    }
};

static ConsoleCls Console;

int pairComparator(const void * a, const void * b)
{
    const Pair pair1 = *(const Pair*)a;
    const Pair pair2 = *(const Pair*)b;
    if ( pair1.dist < pair2.dist ) {
        return -1;
    } else if ( pair1.dist == pair2.dist ) {
        return 0;
    } else {
        return 1;
    }
}

class Heap
{
private:
    Pair* heapArray;
    int maxSize;
public:
    int currentSize;
    Heap(int maxHeapSize = 0)
    {
        maxSize = maxHeapSize;
        currentSize = 0;
        heapArray = new Pair[maxSize];
    }
    bool IsNull() {
        return maxSize == 0;
    }
    bool IsEmpty()
    { return currentSize == 0; }
    bool Insert(Pair value)
    {
        if (currentSize == maxSize)
            return false;

        heapArray[currentSize] = value;
        CascadeUp(currentSize++);
        return true;
    }
    void CascadeUp(int index)
    {
        int parent = (index - 1) / 2;
        Pair bottom = heapArray[index];
        while (index > 0 && heapArray[parent].dist < bottom.dist)
        {
            heapArray[index] = heapArray[parent];
            index = parent;
            parent = (parent - 1) / 2;
        }
        heapArray[index] = bottom;
    }
    Pair Remove() // Remove maximum value Pair
    {
        Pair root = heapArray[0];
        heapArray[0] = heapArray[--currentSize];
        CascadeDown(0);
        return root;
    }
    void CascadeDown(int index)
    {
        int largerChild;
        Pair top = heapArray[index];
        while (index < currentSize / 2)
        {
            int leftChild = 2 * index + 1;
            int rightChild = leftChild + 1;
            if (rightChild < currentSize && heapArray[leftChild].dist < heapArray[rightChild].dist)
                largerChild = rightChild;
            else
                largerChild = leftChild;
            if (top.dist >= heapArray[largerChild].dist)
                break;
            heapArray[index] = heapArray[largerChild];
            index = largerChild;
        }
        heapArray[index] = top;
    }
    bool HeapIncreaseDecreaseKey(int index, Pair newValue)
    {
        if (index < 0 || index >= currentSize)
            return false;
        Pair oldValue = heapArray[index];
        heapArray[index] = newValue;
        if (oldValue.dist < newValue.dist)
            CascadeUp(index);
        else
            CascadeDown(index);
        return true;
    }
    void convertToSorted()
    {
        std::qsort(heapArray, currentSize, sizeof(Pair), pairComparator);
    }
    Pair last()
    {
        return heapArray[currentSize - 1];
    }
    Pair first()
    {
        for (int i = 0; i < currentSize; i++)
            if (heapArray[i].dist != 0)
                return heapArray[i];
        return heapArray[0];
    }
    Pair valueAt(int i)
    {
        return heapArray[i];
    }
    bool Append(Pair value)
    {
        if (currentSize == maxSize)
            return false;

        heapArray[currentSize++] = value;
        return true;
    }


    void DisplayHeap()
    {
        Console.WriteLine();
        Console.Write("Elements of the Heap Array are : ");
        for (int m = 0; m < currentSize; m++)
            if (heapArray[m].dist != 0.0) {
                Console.Write(heapArray[m].dist);
                Console.Write(" ");
            }
            else
                Console.Write("-- ");
        Console.WriteLine();
        int emptyLeaf = 32;
        int itemsPerRow = 1;
        int column = 0;
        int j = 0;
        std::string separator = "...............................";
        Console.WriteLine(separator + separator);
        while (currentSize > 0)
        {
            if (column == 0)
                for (int k = 0; k < emptyLeaf; k++)
                    Console.Write(' ');
            Console.Write(heapArray[j].dist);

            if (++j == currentSize)
                break;
            if (++column == itemsPerRow)
            {
                emptyLeaf /= 2;
                itemsPerRow *= 2;
                column = 0;
                Console.WriteLine();
            }
            else
                for (int k = 0; k < emptyLeaf * 2 - 2; k++)
                    Console.Write(' ');
        }
        Console.WriteLine(std::string("\n") + separator + separator);
    }
};


template<typename T>
T** alloc2d(uint rows, uint cols) {
    T** arr = new T*[rows];
    for (uint i = 0; i < rows; ++i) {
        arr[i] = new T[cols];
    }
    return arr;
}
/*
template<typename T>
T** alloc2d(uint rows, uint cols, const T& default_val) {
    T** arr = new T*[rows];
    for (uint i = 0; i < rows; ++i) {
        arr[i] = new T[cols];
        for (uint j = 0; j < cols; ++j) {
            arr[rowsdefault_val
        }
    }
    return arr;
}*/

class MathCls {
public:
    template<typename T>
    T Sqrt(const T& val) {
        return sqrt(val);
    }

    template<typename T>
    T Floor(const T& val) {
        return floor(val);
    }
};

static MathCls Math;


double partial_seconds(std::chrono::system_clock::duration dur) {
    double milliseconds = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(dur).count());
    return milliseconds / 1000.0;
}


class Program
{
public:
    int cap;
    double ccc;
    int VERBOSE;
    int WW;
    int K;
    Program(int verbose = 1, int ww = 0, int k = 0): VERBOSE(verbose), WW(ww), K(k) {}

    void Main(std::vector<std::string> args)
    {
        double* TS;
        double* xy;
        double* x;
        double* x2;
        double* prevXY;


        double** bsf;
        long** loc1;
        long** loc2;

        double* maX;

        std::chrono::system_clock::time_point Start = std::chrono::system_clock::now();
        std::chrono::system_clock::duration Elapsed;
        int SIZE = std::stoi(args[1]);

        long n = 1;
        long minLength = std::stol(args[2]);
        long maxLength = std::stol(args[3]);
        K = 50;
        if (args.size() > 4)
            K = std::stoi(args[4]);
        ccc = 0.1;
        if (args.size() > 5)
            ccc = std::stod(args[5]);
        WW = SIZE / 2;
        cap = SIZE;// *(int)minLength;

        try
        {
            {
                std::ifstream in(args[0]);
                TS = new double[SIZE];
                x2 = new double[SIZE];
                x = new double[SIZE];
                xy = new double[SIZE];
                prevXY = new double[SIZE];

                bsf = alloc2d<double>(maxLength, K);
                loc1 = alloc2d<long>(maxLength, K);
                loc2 = alloc2d<long>(maxLength, K);

                maX = new double[maxLength];

                TS[0] = 0; x[0] = 0; x2[0] = 0;
                std::string str;
                while ((in >> str) && n < SIZE - 2)
                {
                    double d = std::stod(str);
                    TS[n] = d;
                    x[n] = d + x[n - 1];
                    x2[n] = d * d + x2[n - 1];

                    n++;
                }
            }

            if (VERBOSE == 1) {
                Console.Write("The size of the data is : ");
                Console.WriteLine(n);
            }

            //Initialize bsf
            for (long j = 0; j < maxLength; j++)
                for (long k = 0; k < K; k++)
                    bsf[j][k] = 99999999999999;

            int count = 0;
            int mistakeCount = 0;


            /* The new algorithm for the length independent motif */

            Start = std::chrono::system_clock::now();
            count = 0;

            ////Start: maX computation block

            for (long j = 0; j < maxLength; j++)
            {
                maX[j] = 0;
            }

            for (long i = 1; i < n - maxLength; i++)
            {
                long mLoc = -1;
                long mnLoc = -1;
                double runMax = 0;
                double runMin = 999999999999999.099;


                for (long j = minLength; j < maxLength; j++)
                {

                    double sumY = x[i + j - 1] - x[i - 1];
                    double meanY = sumY / j;
                    double sumY2 = x2[i + j - 1] - x2[i - 1];
                    double sigmaY = Math.Sqrt((sumY2 / j) - (meanY * meanY));




                    if (j == minLength)
                    {
                        runMax = 0;
                        runMin = 999999999999999.099;
                        for (long k = 1; k <= minLength; k++)
                        {
                            double X = (TS[i + k - 1] - meanY) / sigmaY;
                            if (X > runMax)
                            {
                                runMax = X;
                                mLoc = i + k - 1;
                            }
                            if (X < runMin)
                            {
                                runMin = X;
                                mnLoc = i + k - 1;
                            }



                        }

                    }
                    else
                    {
                        runMax = (TS[mLoc] - meanY) / sigmaY;
                        runMin = (TS[mnLoc] - meanY) / sigmaY;
                        double Y = (TS[i + j - 1] - meanY) / sigmaY;
                        if (runMax < Y)
                        {
                            runMax = Y;
                            mLoc = i + j - 1;
                        }
                        if (runMin > Y)
                        {
                            runMin = Y;
                            mnLoc = i + j - 1;
                        }



                    }

                    if (runMax > maX[j])
                        maX[j] = runMax;
                    //if (-runMax*runMin > maX[j])
                    //   maX[j] = -runMax * runMin;
                    if (-runMin > maX[j])
                        maX[j] = -runMin;

                }

                // Console.WriteLine(i+ " " + maX[128]);
            }

            Elapsed = std::chrono::system_clock::now() - Start;

            if (VERBOSE == 1) {
                Console.Write("Time Elapsed till Max Computation : ");
                Console.Write(partial_seconds(Elapsed));
                Console.WriteLine(" ms");
            }

            /// END : maX computation block
            /// 



            Heap List;
            double firstLB = -1;
            for (long j = minLength; j < maxLength; j++)
            {

                //  Console.WriteLine("LENGTH " + j);
                if (List.IsNull() || List.IsEmpty() == true)
                {
                    // List = FindMotifNewSpace(TS, SIZE, j);
                    List = FindMotif(TS, SIZE, j);
                    List.convertToSorted();

                    double dis = List.last().dist;
                    double z = maX[j];
                    double LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                    firstLB = dis / Math.Sqrt(LB);
                    Elapsed = std::chrono::system_clock::now() - Start;

                    int kk = 0;
                    Pair curr = (Pair)List.valueAt(0);
                    loc1[j][kk] = curr.loc1;
                    loc2[j][kk] = curr.loc2;
                    bsf[j][kk] = curr.dist;
                    kk++;
                    for (int nk = 1; nk < List.currentSize && kk < K; nk++)
                    {
                        curr = (Pair)List.valueAt(nk);
                        int sign = 0;
                        for (int mk = 0; mk < kk && kk < K; mk++)
                        {
                            Pair test = (Pair)List.valueAt(mk);
                            // long d1 = j - Math.Abs(curr.loc2 - test.loc2);
                            //long d2 = j - Math.Abs(curr.loc1 - test.loc1);
                            if (isCovering(curr, (int)j, test, (int)j))
                            //if (d1 < 0.2*j || d2 > 0.2*j)
                            {
                                sign = 1;
                                break;
                            }
                        }
                        if (sign == 0)
                        {
                            loc1[j][kk] = curr.loc1;
                            loc2[j][kk] = curr.loc2;
                            bsf[j][kk] = curr.dist;
                            kk++;
                        }
                    }

                    if (VERBOSE == 1) {
                        Console.Write("Time Elapsed till First Motif Computation : ");
                        Console.Write(partial_seconds(Elapsed));
                        Console.WriteLine(" ms");
                    }
                    count++;
                    //printList(List);
                }
                else
                {

                    double LB = 0, best = 0;
                    double dis = List.last().dist;
                    double z = maX[j];
                    LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                    LB = dis / Math.Sqrt(LB);

                    if (VERBOSE == 1)

                        Console.WriteLine(std::to_string(firstLB) + std::string(" FIRSTLB ") + std::to_string(LB) + std::string(" : LB   ") + std::to_string(dis) + std::string(" : dists\n"));
                    Heap *nextList = new Heap(List.currentSize);

                    for (int k = 0; k < List.currentSize; k++)
                    {


                        Pair curr = (Pair)List.valueAt(k);
                        if (curr.loc1 + j >= n || curr.loc2 + j >= n || curr.loc1 == 0 || curr.loc2 == 0)
                            continue;
                        curr.dist = distance(TS, curr.loc1, curr.loc2, j);


                        nextList->Append(curr);


                        // Console.WriteLine(k);
                    }

                    nextList->convertToSorted();

                    int kk = 0;
                    Pair curr1 = (Pair)List.valueAt(0);
                    best = curr1.dist;
                    kk++;
                    for (int nk = 1; nk < List.currentSize && kk < K; nk++)
                    {
                        curr1 = (Pair)List.valueAt(nk);
                        int sign = 0;
                        for (int mk = 0; mk < kk && kk < K; mk++)
                        {
                            Pair test = (Pair)List.valueAt(mk);
                            //  long d1 = j - Math.Abs(curr1.loc2 - test.loc2);
                            // long d2 = j - Math.Abs(curr1.loc1 - test.loc1);
                            if (isCovering(curr1, (int)j, test, (int)j))
                            {
                                sign = 1;
                                break;
                            }
                        }
                        if (sign == 0)
                        {
                            best = curr1.dist;
                            kk++;
                        }
                    }




                    if (best >= firstLB)
                    {

                        List = FindMotif(TS, SIZE, j);
                        List.convertToSorted();
                        dis = List.last().dist;
                        z = maX[j];
                        LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                        firstLB = dis / Math.Sqrt(LB);
                        if (VERBOSE == 1)
                            Console.WriteLine(std::string("\n\nFIRST LB ") + std::to_string(firstLB) + std::string("\n\n"));

                        count++;



                    }
                    else if (best < firstLB)
                    {


                        //Cut down the list to maintain order. || (double)List.GetKey(i) > LB
                        //XXX MEMORY LEAK
                        List = *nextList;
                        z = maX[j];
                        LB = j / (j - 1) + z * z * (j - 1) / (j * j);
                        firstLB = firstLB / Math.Sqrt(LB);

                    }
                    kk = 0;
                    curr1 = (Pair)List.valueAt(0);
                    loc1[j][kk] = curr1.loc1;
                    loc2[j][kk] = curr1.loc2;
                    bsf[j][kk] = curr1.dist;
                    kk++;
                    for (int nk = 1; nk < List.currentSize && kk < K; nk++)
                    {
                        curr1 = (Pair)List.valueAt(nk);
                        int sign = 0;
                        for (int mk = 0; mk < kk && kk < K; mk++)
                        {
                            Pair test = (Pair)List.valueAt(mk);
                            // long d1 = j - Math.Abs(curr1.loc2 - test.loc2);
                            //  long d2 = j - Math.Abs(curr1.loc1 - test.loc1);
                            //  if (d1 > 0.2 * j || d2 > 0.2 * j)
                            if (isCovering(curr1, (int)j, test, (int)j))
                            {
                                sign = 1;
                                break;
                            }
                        }
                        if (sign == 0)
                        {

                            loc1[j][kk] = curr1.loc1;
                            loc2[j][kk] = curr1.loc2;
                            bsf[j][kk] = curr1.dist;
                            kk++;

                        }
                    }

                    Pair P = (Pair)List.valueAt(K - 1);
                    if (VERBOSE == 1)
                    {
                        Console.WriteLine(std::string("Correct Numbers ") + std::to_string(bsf[j][K - 1]) + std::string("\t") + std::to_string(loc1[j][K - 1]) + std::string("\t") + std::to_string(loc2[j][K - 1]));
                        Console.WriteLine(std::to_string(P.dist) + std::string("\t") + std::to_string(P.loc1) + std::string("\t") + std::to_string(P.loc2) + std::string("\t") + std::to_string(j) + std::string("\t") + std::to_string(List.currentSize));
                    }

                }




            }

            Elapsed = std::chrono::system_clock::now() - Start;

            if (VERBOSE == 1)
            {
                Console.Write("Count of FindMotif calls : ");
                Console.WriteLine(count);
                Console.Write("Count of Guesses : ");
                Console.WriteLine(maxLength - minLength - count);
                Console.Write("mistakeCount ");
                Console.WriteLine(mistakeCount);
                Console.Write("Time Elapsed : ");
                Console.Write(partial_seconds(Elapsed));
                Console.WriteLine(" ms");
            }


            int** mark = alloc2d<int>(maxLength, K);

            for (long i = minLength; i < maxLength; i++)
            {

                for (int kk = 0; kk < K; kk++)
                {
                    mark[i][kk] = 0;
                    long j = i;
                    for (int kj = kk + 1; kj < K; kj++)
                        if (isCovering(loc1[i][kk], loc2[i][kk], (int)j, loc1[i][kj], loc2[i][kj], (int)j))
                        {
                            mark[i][kk] = 1;
                            break;
                        }
                    for (j = i + 1; j < maxLength; j++)
                    {
                        for (int kj = 0; kj < K; kj++) {
                            if (isCovering(loc1[j][kj], loc2[j][kj], (int)j, loc1[i][kk], loc2[i][kk], (int)i))
                            {
                                mark[i][kk] = 1;
                                break;

                            }
                        }
                    }
                }
            }

            for (long i = minLength; i < maxLength; i++)
                for (int kk = 0; kk < K; kk++)
                    if (mark[i][kk] == 0) {
                        Console.Write(i);
                        Console.Write(",");
                        Console.Write(loc1[i][kk]);
                        Console.Write(",");
                        Console.Write(loc2[i][kk]);
                        Console.Write(",");
                        Console.Write(bsf[i][kk]);
                        Console.WriteLine();
                    }

        }
        
        catch (std::exception& e)
        {
            Console.WriteLine("The file could not be read:");
            Console.WriteLine(e.what());
        }



    }
/*
    void printList(SortedList L)
    {
        for (int k = 0; k < L.Count; k++)
        {
            Pair PP = (Pair)L.GetByIndex(k);
            Console.WriteLine((double)L.GetKey(k) + "\t" + PP.loc1 + "\t" + PP.loc2);

        }

    }
*/

    Heap FindMotif(double* TS, long SIZE, long LENGTH)
    {

        Heap null_heap;
        try
        {
            double* xy;
            double* x;
            double* x2;
            double* prevXY;


            double bsf;
            long loc1 = -1;
            long loc2 = -1;


            x2 = new double[SIZE];
            x = new double[SIZE];
            xy = new double[SIZE];
            prevXY = new double[SIZE];

            Heap List(cap);

            TS[0] = 0; x[0] = 0; x2[0] = 0;

            long n = 1;

            while (n < SIZE - 2)
            {

                double d = TS[n];
                x[n] = d + x[n - 1];
                x2[n] = d * d + x2[n - 1];
                n++;
            }


            //   Console.WriteLine("The size of the data is : " + n);

            //Initialize bsf
            bsf = 99999999999999;

            // compute the first dot product. Can be done quickly by FFT
            xy[0] = 0;
            for (long i = 1; i <= n - LENGTH; i++)
            {
                xy[i] = 0;
                for (long j = 0; j < LENGTH; j++)
                {
                    xy[i] += TS[i + j] * TS[j + 1];
                }
                prevXY[i] = xy[i];
            }


            // for (long i = 1; i < n - WW - LENGTH; i++)
            for (long i = 1; i < n - LENGTH; i++)
            {
                //process dot product

                if (i > 1)
                {
                    for (long k = i + LENGTH; k <= n - LENGTH; k++)
                    {
                        xy[k] = prevXY[k - 1] - TS[k - 1] * TS[i - 1] + TS[k + LENGTH - 1] * TS[i + LENGTH - 1];
                        prevXY[k - 1] = xy[k - 1];
                    }
                    prevXY[n - LENGTH] = x[n - LENGTH];

                }





                // I have a query at this point which is TS[i:i+j]
                //Find the nearest neighbor
                //for (long k = i + WW + LENGTH; k < n - LENGTH; k++)
                for (long k = i + LENGTH; k < n - LENGTH; k++)
                {

                    long j = LENGTH;
                    //Get the necessary statistics
                    double sumX = x[k + j - 1] - x[k - 1];
                    double meanX = sumX / j;
                    double sumX2 = x2[k + j - 1] - x2[k - 1];
                    double sigmaX = Math.Sqrt((sumX2 / j) - (meanX * meanX));
                    double sumY = x[i + j - 1] - x[i - 1];
                    double meanY = sumY / j;
                    double sumY2 = x2[i + j - 1] - x2[i - 1];
                    double sigmaY = Math.Sqrt((sumY2 / j) - (meanY * meanY));


                    //Compute Distance
                    double corr = (xy[k] - (j * meanX * meanY)) / (j * sigmaX * sigmaY);
                    double dis = Math.Sqrt(2 * j * (1 - (float)corr));

                    //        Console.WriteLine(dis);
                    if (dis < bsf)
                    {
                        bsf = dis;
                        loc1 = i;
                        loc2 = k;
                    }

                    Pair P;
                    P.loc1 = i;
                    P.loc2 = k;
                    P.dist = dis;


                    //Insert in the sorted list
                    if (List.Insert(P) == false)
                    {
                        List.HeapIncreaseDecreaseKey(0, P);

                    }

                }
            }

            //  Console.WriteLine(bsf + "\t" + loc1 + "\t" + loc2 + "\t" + LENGTH);
            return List;


        }
        catch (std::exception& e)
        {
            Console.WriteLine("ERROR in the findMotif function!!!");
            Console.WriteLine("The file could not be read:");
            Console.WriteLine(e.what());
        }
        return null_heap;
    }

/*
    double distanceEAB(double[,] data, long i, long k, long length, double bsf = 999999999.999)
    {
        double sum = 0;
        double bsf2 = bsf * bsf;
        long j;
        for (j = 0; j < length && sum <= bsf2; j++)
        {
            sum += (data[i, j] - data[k, j]) * (data[i, j] - data[k, j]);

        }
        return Math.Sqrt(sum);
    }
*/

    bool isCovering(Pair large, int largeLength, Pair small, int smallLength)
    {
        //assumes loc1 < loc2
        double c = ccc;
        if (large.loc1 - small.loc1 >= 0 && large.loc1 - small.loc1 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc1 + large.loc1 > c * smallLength && largeLength - small.loc1 + large.loc1 <= largeLength)
            return true;
        else if (large.loc2 - small.loc2 >= 0 && large.loc2 - small.loc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc2 + large.loc2 > c * smallLength && largeLength - small.loc2 + large.loc2 <= largeLength)
            return true;
        else if (large.loc1 - small.loc2 >= 0 && large.loc1 - small.loc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc2 + large.loc1 > c * smallLength && largeLength - small.loc2 + large.loc1 <= largeLength)
            return true;
        else if (large.loc1 - small.loc2 >= 0 && large.loc1 - small.loc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - small.loc2 + large.loc1 > c * smallLength && largeLength - small.loc2 + large.loc1 <= largeLength)
            return true;
        else
            return false;
    }
    bool isCovering(long largeLoc1, long largeLoc2, int largeLength, long smallLoc1, long smallLoc2, int smallLength)
    {
        //assumes loc1 < loc2

        double c = ccc;
        if (largeLoc1 - smallLoc1 >= 0 && largeLoc1 - smallLoc1 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc1 + largeLoc1 > c * smallLength && largeLength - smallLoc1 + largeLoc1 <= largeLength)
            return true;
        else if (largeLoc2 - smallLoc2 >= 0 && largeLoc2 - smallLoc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc2 + largeLoc2 > c * smallLength && largeLength - smallLoc2 + largeLoc2 <= largeLength)
            return true;
        else if (largeLoc1 - smallLoc2 >= 0 && largeLoc1 - smallLoc2 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc2 + largeLoc1 > c * smallLength && largeLength - smallLoc2 + largeLoc1 <= largeLength)
            return true;
        else if (largeLoc1 - smallLoc2 >= 0 && largeLoc1 - smallLoc1 < (1 - c) * smallLength)
            return true;
        else if (largeLength - smallLoc2 + largeLoc1 > c * smallLength && largeLength - smallLoc2 + largeLoc1 <= largeLength)
            return true;

        else
            return false;
    }

    double distance(double* TS, long i, long k, long length)
    {
        double xy = 0, x = 0, y = 0, x2 = 0, y2 = 0;
        for (long ii = 0; ii < length; ii++)
        {
            xy += TS[ii + i] * TS[ii + k];
            x += TS[ii + i];
            x2 += TS[ii + i] * TS[ii + i];

            y += TS[ii + k];
            y2 += TS[ii + k] * TS[ii + k];

        }


        double meanX = x / length;
        double sigmaX = Math.Sqrt((x2 / length) - (meanX * meanX));
        double meanY = y / length;
        double sigmaY = Math.Sqrt((y2 / length) - (meanY * meanY));


        //Compute Distance
        double corr = (xy - (length * meanX * meanY)) / (length * sigmaX * sigmaY);
        double dis = Math.Sqrt(2 * length * (1 - (float)corr));
        return dis;

    }

};


int main(int argc, const char* argv[]) {
    Console.WriteLine("Starting");
    Program prog;
    std::vector<std::string>args = {
        std::string("/home/pvnick/Downloads/moen/cpp/LSF5_10.txt"),
        std::string("10000"),
        std::string("128"),
        std::string("256")
    };
    for (auto s: args) {
        Console.WriteLine(s);
    }
    prog.Main(args);
    return 0;
}