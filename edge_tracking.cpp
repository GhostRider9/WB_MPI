#include <stdio.h>
#include <string.h>
#include <math.h>

class EdgeTracker {

    int Width = -1, Height = -1, MaxIter = -1;
    const float cre = -1.36022, cim = 0.0653316, diam = 0.035;
//    float minr = cre - diam * .5, mini = cim - diam * .5;
//    float maxr = cre + diam * .5, maxi = cim + diam * .5;
    float mini, minr, maxr, maxi;
//    const float stepr = (maxr - minr) / Width;
//    const float stepi = (maxi - mini) / Height;
//
    float stepr, stepi;

    enum {
        Loaded = 1, Queued = 2
    };
    int *Data;
    int *Done;
    int QueueSize;
    int *Queue;
    int QueueHead = 0, QueueTail = 0;


    int Iterate(float x, float y) {
        int iter;
        float r = x, i = y;
        for (iter = 0; iter < MaxIter; ++iter) {
            float r2 = r * r, i2 = i * i;
            if (r2 + i2 > 4.0)
                break;
            float ri = r * i;
            i = ri + ri + y;
            r = r2 - i2 + x;
        }
        return iter;
    }

    int inset(float real, float img) {

        // cardioid check
        float img2 = img * img;
        float q = (real - 1.0 / 4.0) * (real - 1.0 / 4.0) + img2;
        if ((q * (q + (real - 1.0 / 4.0))) < (1.0 / 4.0 * img2)) {
            // in cardioid
            return MaxIter;
        }
        // period 2 bulb check
        if (((real + 1) * (real + 1)) + img2 < 1.0 / 16.0) {
            return MaxIter; // in period 2 bulb
        }


        float z_real = real;
        float z_img = img;

        float test_real = z_real;
        float test_img = z_img;
        int period = 8;
        int period_index = 0;
        int iters=0;
        for (iters = 0; iters < MaxIter; iters++) {

            float z2_real = z_real * z_real - z_img * z_img;
            float z2_img = 2.0 * z_real * z_img;
            z_real = z2_real + real;
            z_img = z2_img + img;
            if ((z_real == test_real) && (z_img == test_img)) {
                return MaxIter;
            }
            if (z_real * z_real + z_img * z_img > 4.0) break;
            period_index++;
            if (period_index == period) {
                test_real = z_real;
                test_img = z_img;
                period_index = 0;
                period *= 2;
            }

        }
        return iters;
    }

    void AddQueue(unsigned p) {
        if (Done[p] & Queued) return;
        Done[p] |= Queued;
        Queue[QueueHead++] = p;
        if (QueueHead == QueueSize) QueueHead = 0;
    }

    int Load(unsigned p) {
        if (Done[p] & Loaded) return Data[p];
        unsigned x = p % Width, y = p / Width;
        int result = inset(minr + x * stepr, mini + y * stepi);
        //PutPixel(x, y, result);
        Done[p] |= Loaded;
        return Data[p] = result;
    }


    void Scan(unsigned p) {
        unsigned x = p % Width, y = p / Width;
        int center = Load(p);
        int ll = x >= 1, rr = x < Width - 1;
        int uu = y >= 1, dd = y < Height - 1;
/* ^These are booleans, but bc31 does
* not support the "bool" type */
/* If a neighbor color differs from
* the center, scan the neighbor in turn */
        int l = ll && Load(p - 1) != center;
        int r = rr && Load(p + 1) != center;
        int u = uu && Load(p - Width) != center;
        int d = dd && Load(p + Width) != center;
        if (l) AddQueue(p - 1);
        if (r) AddQueue(p + 1);
        if (u) AddQueue(p - Width);
        if (d) AddQueue(p + Width);
/* The corner pixels (nw,ne,sw,se) are also neighbors */
        if ((uu && ll) && (l || u)) AddQueue(p - Width - 1);
        if ((uu && rr) && (r || u)) AddQueue(p - Width + 1);
        if ((dd && ll) && (l || d)) AddQueue(p + Width - 1);
        if ((dd && rr) && (r || d)) AddQueue(p + Width + 1);
    }
//

public:
    int
    pointsInRegion(float min_real, float max_real, float min_img, float max_img, int maxIterations, int regionWidth,
                   int regionHeight) {


        this->minr = min_real;
        this->maxr = max_real;
        this->mini = min_img;
        this->maxi = max_img;
        MaxIter = maxIterations;
        Width = regionWidth;
        Height = regionHeight;
        stepr = (maxr - minr) / Width;
        printf("maxi: %f mini: %f", this->maxi, this->mini);
        stepi = (maxi - mini) / Height;
        printf("stepr: %f stepi:%f \n", stepr, stepi);
        QueueSize = (Width + Height) * 4;
        this->Queue = new int[QueueSize];

        QueueHead = 0;
        QueueTail = 0;

        this->Data = new int[Width * Height];

        this->Done = new int[Width * Height];

        for (unsigned y = 0; y < Height; y++) {
            AddQueue(y * Width + 0);
            AddQueue(y * Width + (Width - 1));
        }
        for (unsigned x = 1; x < Width - 1; x++) {
            AddQueue(0 * Width + x);
            AddQueue((Height - 1) * Width + x);
        }
        unsigned flag = 0;
        while (QueueTail != QueueHead) {
            unsigned p;
            if (QueueHead <= QueueTail || ++flag & 3) {
                p = Queue[QueueTail++];
                if (QueueTail == QueueSize) QueueTail = 0;
            } else p = Queue[--QueueHead];
            Scan(p);
        }

        for (unsigned p = 0; p < Width * Height; p++) {
            if (Done[p] & Loaded)
                if (!(Done[p + 1] & Loaded)) {
                    // VRAM[p + 1] = VRAM[p],
                    Data[p + 1] = Data[p],
                    Done[p + 1] |= Loaded;
                }
        }
//        printf("computed");
//        for (int i = 0; i < Width; i++) {
//            for (int j = 0; j < Height; j++) {
//                printf("Data: %d\n", Data[Width * j + i]);
//            }
//        }

        int count = 0;
        for (int i = 0; i < Width; i++) {
            for (int j = 0; j < Height; j++) {
                if (Data[j*Width+i] >= MaxIter) {
                    count++;
                }
            }
        }
        return count;


    }

};