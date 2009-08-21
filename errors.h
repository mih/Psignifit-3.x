#ifndef ERRORS_H
#define ERRORS_H

class PsiError {
};

class NotImplementedError : public PsiError {
};

class BadArgumentError : public PsiError {
};

class BadIndexError : public PsiError {
};

#endif
