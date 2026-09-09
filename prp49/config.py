"""Configuration loading (YAML with nested access)."""
import os
import yaml


class Config:
    def __init__(self, d):
        self._d = d

    @classmethod
    def from_yaml(cls, path):
        with open(path, encoding='utf-8') as f:
            d = yaml.safe_load(f)
        return cls(d)

    def __getattr__(self, name):
        v = self._d.get(name)
        if isinstance(v, dict):
            return Config(v)
        return v

    def get(self, name, default=None):
        return self._d.get(name, default)

    def resolve(self, path, base_dir=None):
        """Resolve relative path against the config file's directory."""
        if not path:
            return path
        if os.path.isabs(path):
            return path
        if base_dir is None:
            base_dir = getattr(self, '_base_dir', '.')
        return os.path.normpath(os.path.join(base_dir, path))
