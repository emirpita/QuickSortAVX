void scalar_partition_epi32(uint32_t* array, const uint32_t pivot, int& left, int& right) {

    while (left <= right) {

        while (array[left] < pivot) {
            left += 1;
        }

        while (array[right] > pivot) {
            right -= 1;
        }

        if (left <= right) {
            const uint32_t t = array[left];
            array[left]      = array[right];
            array[right]     = t;

            left  += 1;
            right -= 1;
        }
    }
}
void scalar_partition_ps(float* array, const float pivot, int& left, int& right) {

    while (left <= right) {

        while (array[left] < pivot) {
            left += 1;
        }

        while (array[right] > pivot) {
            right -= 1;
        }

        if (left <= right) {
            const float t = array[left];
            array[left] = array[right];
            array[right] = t;

            left += 1;
            right -= 1;
        }
    }
}

void scalar_partition_16(uint16_t *array, const uint16_t pivot, int &left, int &right) {

    while (left <= right) {

        while (array[left] < pivot) {
            left += 1;
        }

        while (array[right] > pivot) {
            right -= 1;
        }

        if (left <= right) {
            const uint16_t t = array[left];
            array[left] = array[right];
            array[right] = t;

            left += 1;
            right -= 1;
        }
    }
}

void scalar_partition_64(uint64_t *array, const uint64_t pivot, int &left, int &right) {

    while (left <= right) {

        while (array[left] < pivot) {
            left += 1;
        }

        while (array[right] > pivot) {
            right -= 1;
        }

        if (left <= right) {
            const uint64_t t = array[left];
            array[left] = array[right];
            array[right] = t;

            left += 1;
            right -= 1;
        }
    }
}
void scalar_partition_8(uint8_t *array, const uint8_t pivot, int &left, int &right) {

    while (left <= right) {

        while (array[left] < pivot) {
            left += 1;
        }

        while (array[right] > pivot) {
            right -= 1;
        }

        if (left <= right) {
            const uint8_t t = array[left];
            array[left] = array[right];
            array[right] = t;

            left += 1;
            right -= 1;
        }
    }
}


void scalar_partition_pd(double *array, const double pivot, int &left, int &right) {

    while (left <= right) {

        while (array[left] < pivot) {
            left += 1;
        }

        while (array[right] > pivot) {
            right -= 1;
        }

        if (left <= right) {
            const double t = array[left];
            array[left] = array[right];
            array[right] = t;

            left += 1;
            right -= 1;
        }
    }
}


