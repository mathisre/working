#ifndef treeutilsH
#define treeutilsH


class devtree
{
 public:
  int val;
  static devtree *hasdev(devtree *root, int val);
  static devtree *remove(devtree *dt);
  static devtree *add(devtree *dt, int val);
  devtree(int setval) { less=NULL; more=NULL; parent = NULL; val=setval;}
  static bool equal(devtree *root1, devtree *root2);
  static void deltree(devtree* root);
  static int *vals(devtree *root, int size);
  static devtree *successor(devtree *start);

  // private:
  static devtree *findmin(devtree* root);
  devtree *less;
  devtree *more;
  devtree *parent;
};

#endif
