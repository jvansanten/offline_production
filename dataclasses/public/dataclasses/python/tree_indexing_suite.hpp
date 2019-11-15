//  tree_indexing_suite: Python-land indexing and manipulation for new-style I3Tree

#ifndef DATACLASSES_PYTHON_TREE_INDEXING_SUITE_HPP_INCLUDED
#define DATACLASSES_PYTHON_TREE_INDEXING_SUITE_HPP_INCLUDED

#include <boost/python/iterator.hpp>
#include <boost/type_index.hpp>
#include <icetray/python/get_class.hpp>
#include <dataclasses/physics/I3MCTreeUtils.h>

namespace TreeBase {

template <class Container>
class tree_indexing_suite : public boost::python::def_visitor<tree_indexing_suite<Container > >{
public:

    typedef typename Container::value_type value_type;
    typedef typename Container::key_type key_type;
    typedef typename Container::optional_value optional_value;
    typedef typename Container::pre_order_iterator pre_order_iterator;
    typedef typename Container::post_order_iterator post_order_iterator;
    typedef typename Container::sibling_iterator sibling_iterator;

    struct not_found_exception : std::exception
    {
      std::string msg;
      not_found_exception(std::string m) : msg(m) { }
      ~not_found_exception() throw () { }
      char const* what() const  throw() { return msg.c_str(); }
    };

    static void translate(const not_found_exception& e)
    {
        // Use the Python 'C' API to set up an exception object
        PyErr_SetString(PyExc_IndexError, e.what());
    }

    static value_type getItem(const optional_value& ptr, std::string err="")
    {
      if (ptr)
        return *ptr;
      else
        throw not_found_exception(err);
    }

    static value_type get_best(const Container& tree, boost::python::object func) {
      return getItem(I3MCTreeUtils::GetBest(tree, func),"no items in tree");
    }

    static const std::vector<value_type> get_filter(const Container& tree, boost::python::object func) {
      return I3MCTreeUtils::GetFilter(tree, func);
    }

    static value_type get_best_filter(const Container& tree, boost::python::object func1,boost::python::object func2) {
      return getItem(I3MCTreeUtils::GetBestFilter(tree, func1, func2),"no items found");
    }

    static value_type get_head(const Container& t)
    { return getItem(t.get_head(),"no head in tree"); }

    static value_type parent(const Container& t,const key_type& p)
    { return getItem(t.parent(p),"item not found or no parent"); }

    static value_type previous_sibling(const Container& t,const key_type& p)
    { return getItem(t.previous_sibling(p),"item not found or no sibling"); }

    static value_type next_sibling(const Container& t,const key_type& p)
    { return getItem(t.next_sibling(p),"item not found or no sibling"); }

    static value_type first_child(const Container& t,const key_type& p)
    { return getItem(t.first_child(p),"item not found or no chilren"); }

    static value_type* at(Container& t,const key_type& p)
    {
      typename Container::iterator iter(t,p);
      if (iter == t.end())
        throw not_found_exception("item not found in tree");
      else
        return &(*iter);
    }

    static value_type* at2(Container& t,int num)
    {
      value_type* ptr = NULL;
      if (num < 0) {
        if (unsigned(-num) > t.size())
          throw not_found_exception("item index not found in tree");
        num = num+t.size();
      } else {
        if (unsigned(num) >= t.size())
          throw not_found_exception("item index not found in tree");
      }
      typename Container::iterator iter = t.begin();
      for(int cnt=0;iter != t.end();cnt++) {
        if (cnt == num) {
          ptr = &(*iter);
          break;
        }
        iter++;
      }
      if (ptr != NULL)
        return ptr;
      else
        throw not_found_exception("item index not found in tree");
    }

    static std::vector<value_type> at3(Container& t,boost::python::slice sl)
    {
      std::vector<value_type> ret;
      int start = 0,stop=t.size(),step=1;
      if (sl.start() != boost::python::object())
        start = boost::python::extract<int>(sl.start());
      if (sl.stop() != boost::python::object()) {
        stop = boost::python::extract<int>(sl.stop());
        if (stop < 0)
            stop += t.size();
      }
      if (sl.step() != boost::python::object())
        step = boost::python::extract<int>(sl.step());
      typename Container::iterator iter = t.begin();
      for(int cnt=0;iter != t.end();cnt++) {
        if (cnt >= start && cnt < stop && cnt % step == 0)
          ret.push_back(*iter);
        if (cnt >= stop)
          break;
        iter++;
      }
      if (ret.size() > 0 && step < 1)
        std::reverse(ret.begin(),ret.end());
      return ret;
    }

    struct outer {
      static inline boost::python::object pass_through(boost::python::object const& o) { return o; }
      class sib_iter {
        public:
          sib_iter(const Container& t, const key_type& p)
            : i_(t,p), end_(t) { end_ = t.cend_sibling(); }
          value_type next()
          {
            if (i_ == end_) {
              PyErr_SetString(PyExc_StopIteration, "No more data.");
              boost::python::throw_error_already_set();
            }
            return *i_++;
          }
        private:
          typename Container::sibling_const_iterator i_;
          typename Container::sibling_const_iterator end_;
      };
    };

    static typename outer::sib_iter sibling_iter(const Container& t, const key_type& p)
    { return typename outer::sib_iter(t,p); }

    static bool nonempty(const Container &t)
    { return !t.empty(); }

    static bool contains(const Container& t,const key_type& p)
    { return bool(t.at(p)); }

    template <typename T>
    static std::string type_name() { return boost::typeindex::type_id<T>().pretty_name(); }

    template <class Class>
    static void
    visit(Class& cl)
    {
        boost::python::register_exception_translator<not_found_exception>(&translate);

        boost::python::scope outer = cl;

        const std::vector<value_type> (Container::*children)(const key_type&) const  = &Container::children;
        void (Container::*erase)(const key_type&) = &Container::erase;
        void (Container::*erase_children)(const key_type&) = &Container::erase_children;
        void (Container::*append_child_subtree)(const key_type&,const TreeBase::Tree<value_type,key_type>&,const key_type&) = &Container::append_child;
        void (Container::*append_children)(const key_type&,const std::vector<value_type>&) = &Container::append_children;
        void (Container::*insert_head)(const value_type&) = &Container::insert;
        void (Container::*insert_after_head)(const value_type&) = &Container::insert_after;
        void (Container::*insert)(const key_type&,const value_type&) = &Container::insert;
        void (Container::*insert_after)(const key_type&,const value_type&) = &Container::insert_after;
        void (Container::*insert_subtree)(const key_type&,const TreeBase::Tree<value_type,key_type>&,const key_type&) = &Container::insert_subtree;
        void (Container::*insert_subtree_after)(const key_type&,const TreeBase::Tree<value_type,key_type>&,const key_type&) = &Container::insert_subtree_after;
        void (Container::*replace)(const key_type&,const value_type&) = &Container::replace;
        void (Container::*replace_subtree)(const key_type&,const TreeBase::Tree<value_type,key_type>&,const key_type&) = &Container::replace;
        void (Container::*flatten)(const key_type&) = &Container::flatten;
        void (Container::*reparent)(const key_type&,const key_type&) = &Container::reparent;
        typename Container::size_type (Container::*depth)(const key_type&) const = &Container::depth;
        typename Container::size_type (Container::*number_of_children)(const key_type&) const = &Container::number_of_children;
        typename Container::size_type (Container::*number_of_siblings)(const key_type&) const = &Container::number_of_siblings;
        bool (Container::*is_in_subtree)(const key_type&,const key_type&) const = &Container::is_in_subtree;
        bool (Container::*subtree_in_tree)(const TreeBase::Tree<value_type,key_type>&,const key_type&) const = &Container::subtree_in_tree;

        cl
            // I3MCTreeUtils
            .def("get_best", get_best, ("Get the best matching "+type_name<value_type>()).c_str())
            .def("get_filter", &get_filter, ("Get the "+type_name<value_type>()+"s passing the filter").c_str())
            .def("get_best_filter", &get_best_filter, ("Get the best matching "+type_name<value_type>()+" passing the filter").c_str())

            // Base Class Methods
            .def("get_head", &get_head, "Get the left-most primary (the root or head of the tree)")
            .def("get_heads", &Container::get_heads, "Get a list of all primaries (the roots or heads of the tree)")
            .def("at", &at, boost::python::return_internal_reference<>(), ("Get the I3Particle represented by the "+type_name<key_type>()).c_str())
            .def("parent", &parent, ("Get the parent of the "+type_name<key_type>()).c_str())
            .def("previous_sibling", &previous_sibling, ("Get the previous sibling of the "+type_name<key_type>()).c_str())
            .def("next_sibling", &next_sibling, ("Get the next sibling of the "+type_name<key_type>()).c_str())
            .def("children", children, ("Get the children of the "+type_name<key_type>()).c_str())
            .def("first_child", &first_child, ("Get the first (left-most) child of the "+type_name<key_type>()).c_str())
            .def("clear", &Container::clear, "Clear everything from the tree")
            .def("erase", erase, ("Erase the "+type_name<key_type>()+" and all children").c_str())
            .def("erase_children", erase_children, ("Erase only the children of the "+type_name<key_type>()+" (keeping the "+type_name<key_type>()+" itself)").c_str())
            .def("append_child", append_child_subtree)
            .def("append_children", append_children, ("Add multiple children to an "+type_name<key_type>()).c_str())
            .def("insert", insert_head, ("Add an "+type_name<value_type>()+" at the root level, before other "+type_name<key_type>()+"s").c_str())
            .def("insert_after", insert_after_head, ("Add an "+type_name<value_type>()+" at the root level, after other "+type_name<key_type>()+"s").c_str())
            .def("insert", insert, ("Add an "+type_name<value_type>()+" before the sibling "+type_name<key_type>()).c_str())
            .def("insert_after", insert_after, ("Add an "+type_name<value_type>()+" after the sibling "+type_name<key_type>()).c_str())
            .def("insert_subtree", insert_subtree, ("Add a subtree of "+type_name<value_type>()+" before the sibling "+type_name<key_type>()).c_str())
            .def("insert_subtree_after", insert_subtree_after, ("Add a subtree of "+type_name<value_type>()+" after the sibling "+type_name<key_type>()).c_str())
            .def("replace", replace, ("Replace an "+type_name<key_type>()+" with another "+type_name<value_type>()).c_str())
            .def("replace", replace_subtree, ("Replace an "+type_name<key_type>()+" and all children with another subtree of "+type_name<value_type>()+"s").c_str())
            .def("flatten", flatten, ("Move the children of the "+type_name<value_type>()+" to be siblings after it").c_str())
            .def("reparent", reparent, ("Move all children to an "+type_name<key_type>()+", from another "+type_name<key_type>()+" in the tree").c_str())
            .def("merge",&Container::merge, "Merge two trees, modifying the first tree")
            .def("size", &Container::size, ("Get the number of "+type_name<value_type>()+"s in the tree").c_str())
            .def("empty", &Container::empty, "Is the tree empty?")
            .def("swap", &Container::swap, "Swap the contents of another tree with this one.")
            .def("depth", depth, ("Get the depth from the "+type_name<key_type>()+" to the primary").c_str())
            .def("number_of_children", number_of_children, ("Get the number of children an "+type_name<key_type>()+" has").c_str())
            .def("number_of_siblings", number_of_siblings, ("Get the number of siblings an "+type_name<key_type>()+" has").c_str())
            .def("is_in_subtree", is_in_subtree, ("Is an "+type_name<value_type>()+" in a subtree?").c_str())
            .def("subtree_in_tree", subtree_in_tree, "Is any part of a subtree in the tree?")

            // Iterators
            .def("pre_order_iter", boost::python::range<boost::python::return_value_policy<boost::python::copy_const_reference> >
              (
                (typename Container::pre_order_const_iterator(Container::*)() const) &Container::cbegin, 
                (typename Container::pre_order_const_iterator(Container::*)() const) &Container::cend
              ),
              "Pre order iterator for tree. This is the default iterator."
            )
            .def("post_order_iter", boost::python::range<boost::python::return_value_policy<boost::python::copy_const_reference> >
              (
                (typename Container::post_order_const_iterator(Container::*)() const) &Container::cbegin_post, 
                (typename Container::post_order_const_iterator(Container::*)() const) &Container::cend_post
              ),
              "Post order iterator for tree"
            )
            .def("sibling_iter",sibling_iter, "Sibling iterator for tree. Takes a key, provides all next_siblings of it.")
            .def("fast_iter", boost::python::range<boost::python::return_value_policy<boost::python::copy_const_reference> >
              (
                (typename Container::fast_const_iterator(Container::*)() const) &Container::cbegin_fast, 
                (typename Container::fast_const_iterator(Container::*)() const) &Container::cend_fast
              ),
              "Fast iterator for tree. Fast but unordered traversal."
            )
            .def("leaf_iter", boost::python::range<boost::python::return_value_policy<boost::python::copy_const_reference> >
              (
                (typename Container::leaf_const_iterator(Container::*)() const) &Container::cbegin_leaf, 
                (typename Container::leaf_const_iterator(Container::*)() const) &Container::cend_leaf
              ),
              "Leaf iterator for tree. Fast but unordered leaf traversal."
            )

            // Python Special Methods
            .def("__nonzero__", &nonempty)
            .def("__bool__", &nonempty)
            .def("__len__", &Container::size)
            .def("__iter__", boost::python::range<boost::python::return_value_policy<boost::python::copy_const_reference> >
              (
                (typename Container::const_iterator(Container::*)() const) &Container::cbegin, 
                (typename Container::const_iterator(Container::*)() const) &Container::cend
              )
            )
            .def("__contains__", &contains)
            .def("__getitem__", &at, boost::python::return_internal_reference<>())
            .def("__getitem__", &at2, boost::python::return_internal_reference<>())
            .def("__getitem__", &at3)
            .def("__delitem__", erase)

            .def("__value_type__", boost::python::get_class<I3MCTree::value_type>)
            .staticmethod("__value_type__")
        ;

        boost::python::class_<typename outer::sib_iter>("_Sibling_Iter_",
          "Sibling iterator object for tree. DO NOT CALL DIRECTLY. Instead use the sibling_iter method",
          boost::python::no_init)
          .def("next",&outer::sib_iter::next)
          .def("__next__",&outer::sib_iter::next)
          .def("__iter__",&outer::pass_through)
        ;
    }
    
};


} // namespace TreeBase

#endif
