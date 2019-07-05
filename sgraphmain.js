// SGraph GUI support
    // Prototypes
    if (!String.prototype.startsWith) {
        String.prototype.startsWith = function (searchString, position) {
            return this.substr(position || 0, searchString.length) === searchString;
        };
    }
    if (!String.prototype.endsWith) {
        String.prototype.endsWith = function (search, this_len) {
            if (this_len === undefined || this_len > this.length) {
                this_len = this.length;
            }
            return this.substring(this_len - search.length, this_len) === search;
        };
    }
    var sg_main = sgraphgui();
    
    function BlockMove(event) {
        // Tell Safari not to move the window.
        event.preventDefault();
    }
