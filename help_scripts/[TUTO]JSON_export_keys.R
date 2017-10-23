swagger <- fromJSON('
  {
                    "Key 1 at level 1": "value x",
                    "Key 2 at level 1": "value x",
                    "Key 3 at level 1": {
                    "Key 1 at level 2": {
                    "Key 1 at level 3": "value x",
                    "Key 2 at level 3": "value x",
                    "Key 3 at level 3": "value x"
                    },
                    "Key 2 at level 2": {
                    "Key 4 at level 3": "value x",
                    "Key 5 at level 3": "value x",
                    "Key 6 at level 3": "value x"
                    }
                    }
                    }
                    ')

map(swagger, names) %>% unlist

all_keys <- map(swagger, names)
all_keys
all_keys$"Key 3 at level 1"
