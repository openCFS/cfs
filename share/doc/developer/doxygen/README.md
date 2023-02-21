doxygen {#doxygen_README}
=========================

Usage Tipps
----------------
* read the [doxygen manual](https://doxygen.nl/manual)
* look at the examples given in [main.md](@ref mainpage)


Bootstrap Design
----------------

### General approach
we use the doxygen generated html and put a bootstrap-design on top.
This works by
1. setting a custom header and footer to include
 * [bootstrap](https://getbootstrap.com) js/css
 * external [jQuery](https://jquery.com/) (because the one used by doxygen is too old for recent bootstrap)
2. applying bootrstrap styling via javascript in `doxy-boot.js`
3. our own css `customdoxygen.css` (which shouls be used as sparesly as possible)

### External Dependencies and Updates
Bootstrap is currently at version `3.3.7`, and jQuery at version `3.2.1`.
For development is makes sense to include css/js form online sources:
 * https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css
 * https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js
 * https://code.jquery.com/jquery-3.2.1.min.js
This way the versions are clear and one can rapidly test different versions.
For production use on our webserver we use a local copy.

### Doxygen Output
Doxygen output is controlled by the config file `doxy-config.in`.
For all config options see the [documentation](https://doxygen.nl/manual/).

The output can be influenced by a `DoxygenLayout.xml` file.
The default is created by `doxygen -l`.
It can be modified and used in the dogycen config (`doxy-config.in`).
 
#### Class page
The class page has the following elements:
- `<div id="main-nav">` smartmenu syle menu (Main Page, Related Pages, Class, ...)
- `<div class="header">` navigation quicklinks in page
- `<div class="contents">` class description


### Styling of Specific Doxygen Elements

#### Navbar and Menu
One can add custom items via doxygen's configurable [layout](https://doxygen.nl/manual/customize.html#layout).
This allows to add specific menur items.

The [bootstrap navbar](http://getbootstrap.com/docs/3.3/components/#navbar) is then defined in `header.html.in`.
The menu items are produced by doxygen (>=1.1.13) by javascript [jQuery smartmenus](https://www.smartmenus.org/).
It can be directly incorporated into a bootstrap design by the [bootstrap-addon](http://vadikom.github.io/smartmenus/src/demo/bootstrap-navbar.html).
Some modifications were necessary by our custom `doxy-boot.js`

#### Bootstrap Grid
It seems to be difficult to use the grid system e.g. in class definitions, because doxygen does not put everything in fitting `<div>`, e.g. `<h2>` followed by `<div class="textblock">` should be together.
This seems very tedious to collect via js.
Therefore, we do not use anything else than a 12-column container.


### References & Further Information
* A nice example is [biogears](https://biogearsengine.com/documentation/index.html).
* Idea of `doxy-boot.js` from [doxy-bootstrap](https://github.com/StratifyLabs/Doxygen-Bootstrap)
* Part of the concept is from [stratifylabs.co](https://stratifylabs.co/embedded%20design%20tips/2014/01/07/Tips-Integrating-Doxygen-and-Bootstrap/)
* for jQuery functions used in `doxy-boot.js` see the [jQuery-docs](https://api.jquery.com/category/manipulation/)
* bootstrap [styling](https://getbootstrap.com/docs/3.3/css/), [components](https://getbootstrap.com/docs/3.3/components/) or [examples](https://getbootstrap.com/docs/3.3/getting-started/#examples).
* [bootstrap grid system](https://getbootstrap.com/docs/3.3/css/#grid)


### ToDo
* ~~table in function description: class fieldtable~~
* ~~format class=memberdecls~~
* ~~fix toggleInherit~~
* ~~format doxygen code~~
* ~~reference list format~~
* ~~toggleLevel~~
* ~~class index table format~~
* ~~class list~~
* ~~empty row in `classCoupledField_1_1ADBInt.html`~~
* ~~single class view header colums: hadling left, info small and right~~
* ~~format 'improvements'~~
* ~~make bullet points appear again in lists (in `div.textblock` = markdown)~~
* ~~ include bootstrap ~~
* define menu items: [page order](https://stackoverflow.com/questions/18001897/how-do-you-order-doxygen-custom-pages)
1. menu second level items in mobile view are not correct
4. add SEARCH by doxygen
8. cleanup the files `doc-doxygen.cmake.in` and `share/doc/developer/CMakeLists.txt`
