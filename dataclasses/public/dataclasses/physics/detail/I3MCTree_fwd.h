#ifndef I3MCTREE_FWD_H
#define I3MCTREE_FWD_H

namespace TreeBase {
  template<typename T> class TreeNode;
  template<typename T, typename Key=T, typename Hash=hash<Key> > class Tree;
}

#endif
