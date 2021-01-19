#include <limits>

template<typename NumericType>
class InputData {

protected:
    NumericType *array;
    size_t n;

public:
    InputData(size_t count) : n(count) {
        assert(n > 0);
        array = new NumericType[n];
    }

    ~InputData() {
        delete[] array;
    }

public:
    NumericType *pointer() {
        return array;
    }

    size_t count() const {
        return n;
    }

    size_t size() const {
        return n * sizeof(NumericType);
    }
};

template<typename NumericType>
class InputAscending : public InputData<NumericType> {

    using super = InputData<NumericType>;

public:
    InputAscending(size_t count) : super(count) {
        for (size_t i = 0; i < InputData<NumericType>::n; i++) {
            InputData<NumericType>::array[i] = NumericType(i % std::numeric_limits<NumericType>::max());
        }
    }
};

template<typename NumericType>
class InputDescending : public InputData<NumericType> {

    using super = InputData<NumericType>;

public:
    InputDescending(size_t count) : super(count) {
        for (size_t i = 0; i < InputData<NumericType>::n; i++) {
            InputData<NumericType>::array[i] = NumericType(fmod(
                    (InputData<NumericType>::n - i + 1) ,std::numeric_limits<NumericType>::max()));
        }
    }
};

template<typename NumericType>
class InputRandomFew : public InputData<NumericType> {

    using super = InputData<NumericType>;

public:
    InputRandomFew(size_t count) : super(count) {
        for (size_t i = 0; i < InputData<NumericType>::n; i++) {
            InputData<NumericType>::array[i] = NumericType((rand() % 10) % std::numeric_limits<NumericType>::max());
        }
    }
};

template<typename NumericType>
class InputRandom : public InputData<NumericType> {

    using super = InputData<NumericType>;

public:
    InputRandom(size_t count) : super(count) {
        for (size_t i = 0; i < InputData<NumericType>::n; i++) {
            InputData<NumericType>::array[i] = NumericType((rand()) % std::numeric_limits<NumericType>::max());
        }
    }
};

template<typename NumericType>
class InputRandomUnique : public InputData<NumericType> {

    using super = InputData<NumericType>;

public:
    InputRandomUnique(size_t count) : super(count) {
        for (size_t i = 0; i < InputData<NumericType>::n; i++) {
            InputData<NumericType>::array[i] = NumericType(i % std::numeric_limits<NumericType>::max());
        }

        shuffle();
    }

private:
    void shuffle() {
        for (size_t i = 0; i < InputData<NumericType>::n; i++) {
            size_t j = rand() % (InputData<NumericType>::n - i);

            const NumericType t = InputData<NumericType>::array[i];
            InputData<NumericType>::array[i] = InputData<NumericType>::array[j];
            InputData<NumericType>::array[j] = t;
        }
    }
};
